"""Production task catalog, parser-backed forms and safe subprocess commands."""

from __future__ import annotations

import argparse
from contextlib import redirect_stderr
from dataclasses import dataclass
from functools import lru_cache
import importlib.util
from io import StringIO
import json
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Callable, TypedDict, cast

SRC = Path(__file__).resolve().parents[1]
Value = str | bool


@dataclass(frozen=True)
class Task:
    name: str
    step: str
    title: str

    @property
    def script(self) -> Path:
        return SRC / self.step / f"{self.name}.py"


TASKS: tuple[Task, ...] = (
    Task("prepare_states", "0001-cyclic-peptide-preparation", "准备环肽化学状态"),
    Task("generate_conformers", "0001-cyclic-peptide-preparation", "生成并筛选环肽构象"),
    Task("select_receptor", "0002-receptor-gprotein-preparation", "选择受体复合物"),
    Task("repair_receptor", "0002-receptor-gprotein-preparation", "补建缺失结构"),
    Task("receptor_qc", "0002-receptor-gprotein-preparation", "检查受体几何"),
    Task("refine_receptor_tails", "0002-receptor-gprotein-preparation", "调整未解析尾段"),
    Task("prepare_receptor", "0002-receptor-gprotein-preparation", "准备微状态与显式氢"),
    Task("orient_receptor", "0002-receptor-gprotein-preparation", "整体定向"),
    Task("receptor_to_pdbqt", "0002-receptor-gprotein-preparation", "导出去胆固醇受体 PDBQT"),
    Task("dock_ligand", "0003-binding-site-selection", "单构象 Vina 对接"),
    Task("run_pipeline", "0003-binding-site-selection", "批量对接与重原子复合物筛选"),
)
STEP_TITLES: dict[str, str] = {
    "0001-cyclic-peptide-preparation": "0001  环肽准备",
    "0002-receptor-gprotein-preparation": "0002  受体与 G 蛋白准备",
    "0003-binding-site-selection": "0003  结合位点筛选",
}
LABELS: dict[str, str] = {
    "input": "输入结构", "input_dir": "输入目录", "output": "输出文件",
    "output_dir": "输出目录", "profile": "体系配置 JSON", "states": "化学状态",
    "template": "观测模板 PDB", "fasta": "完整序列 FASTA", "model": "待检查模型",
    "transform": "坐标变换 JSON", "source_cif": "来源 mmCIF", "forcefield": "力场 XML",
    "toppar": "CHARMM toppar 归档", "chemistry_dir": "配套化学文件目录",
    "receptor": "受体 PDBQT", "ligand": "环肽 SDF", "center": "搜索盒中心 / Å",
    "size": "搜索盒全边长 / Å", "record": "SDF 记录索引（从 0 起）",
    "seed": "随机种子", "threads": "嵌入线程数", "cpu_threads": "CPU 线程数",
    "cpu": "Vina CPU 数", "ph": "pH", "models": "候选模型数", "loops": "每模型 loop 数",
    "attempts_per_state": "每状态采样数", "keep": "代表数上限", "seed_batches": "种子批次数",
    "exhaustiveness": "搜索投入", "max_evals": "每轨迹评估预算", "n_poses": "姿态数上限",
    "inverse": "使用逆变换",
    "config": "对接流程配置 JSON", "resume": "核验后续跑", "screen_only": "仅筛选已有对接结果",
    "limit_jobs": "任务数上限（0 为全部）", "plan_only": "仅检查输入与任务计划",
}


def get_task(name: str) -> Task:
    for task in TASKS:
        if task.name == name:
            return task
    raise ValueError(f"未知任务：{name}")


@lru_cache(maxsize=None)
def get_parser(name: str) -> argparse.ArgumentParser:
    task = get_task(name)
    spec = importlib.util.spec_from_file_location(f"_clemain_{name}", task.script)
    if spec is None or spec.loader is None:
        raise ValueError(f"无法读取入口：{task.script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    builder = cast(Callable[[], argparse.ArgumentParser], getattr(module, "build_parser"))
    parser = builder()
    parser.prog = task.script.name
    return parser


@dataclass(frozen=True)
class Field:
    action: argparse.Action

    @property
    def option(self) -> str:
        return next(option for option in self.action.option_strings if option.startswith("--"))

    @property
    def label(self) -> str:
        return LABELS.get(self.action.dest, self.option)

    @property
    def flag(self) -> bool:
        return self.action.nargs == 0

    @property
    def initial(self) -> Value:
        value = self.action.default
        if self.flag:
            return bool(value)
        if value is None:
            return ""
        if isinstance(value, (tuple, list)):
            return shlex.join(str(item) for item in value)
        return str(value)

    @property
    def hint(self) -> str:
        if self.flag:
            return "Enter 或空格切换开关"
        if self.action.type is Path:
            return "路径可含空格和中文；相对路径以工作目录为准，支持 ~"
        if self.action.nargs is not None:
            count = f"{self.action.nargs} 个" if isinstance(self.action.nargs, int) else "一个或多个"
            return f"输入{count}值，以空格分隔；{self.action.help or ''}"
        return self.action.help or "留空的可选项使用生产脚本默认值"


def fields_for(name: str) -> list[Field]:
    return [Field(action) for action in get_parser(name)._actions
            if action.option_strings and "--help" not in action.option_strings]


@dataclass
class Job:
    task: str
    workdir: Path
    values: dict[str, Value]


def new_job(name: str, workdir: Path) -> Job:
    return Job(name, workdir.expanduser().resolve(), {field.option: field.initial for field in fields_for(name)})


def build_command(job: Job, *, check_paths: bool = True) -> list[str]:
    task = get_task(job.task)
    fields = fields_for(job.task)
    known = {field.option for field in fields}
    if unknown := set(job.values) - known:
        raise ValueError(f"未知参数：{', '.join(sorted(unknown))}")
    args: list[str] = []
    for field in fields:
        value = job.values.get(field.option, field.initial)
        if field.flag:
            if not isinstance(value, bool):
                raise ValueError(f"{field.option} 必须是布尔值")
            if value:
                args.append(field.option)
        else:
            if not isinstance(value, str):
                raise ValueError(f"{field.option} 必须是文本")
            if not value.strip():
                continue
            if field.action.nargs is not None:
                args.extend([field.option, *shlex.split(value)])
            else:
                text = str(Path(value).expanduser()) if field.action.type is Path else value
                args.append(f"{field.option}={text}")
    errors = StringIO()
    try:
        with redirect_stderr(errors):
            namespace = get_parser(job.task).parse_args(args)
    except SystemExit as error:
        message = errors.getvalue().strip().splitlines()
        raise ValueError(message[-1] if message else "参数解析失败") from error
    if check_paths:
        if not job.workdir.is_dir():
            raise ValueError(f"工作目录不存在：{job.workdir}")
        for field in fields:
            value = getattr(namespace, field.action.dest)
            if not isinstance(value, Path):
                continue
            path = value if value.is_absolute() else job.workdir / value
            if field.option in {"--output", "--output-dir"}:
                if path.exists() or path.is_symlink():
                    pipeline_preview = job.task == "run_pipeline" and namespace.plan_only
                    pipeline_resume = (job.task == "run_pipeline" and namespace.resume and path.is_dir()
                                       and (path / "run_manifest.json").is_file())
                    if not (pipeline_preview or pipeline_resume):
                        raise ValueError(f"输出已存在，请使用新路径：{path}")
            elif field.option in {"--input-dir", "--chemistry-dir"}:
                if not path.is_dir():
                    raise ValueError(f"{field.option} 目录不存在：{path}")
            elif not path.is_file():
                raise ValueError(f"{field.option} 文件不存在：{path}")
    return [sys.executable, "-u", str(task.script), *args]


def command_text(job: Job, command: list[str]) -> str:
    return f"cd -- {shlex.quote(str(job.workdir))}\n{shlex.join(command)}"


class SavedJob(TypedDict):
    version: int
    task: str
    workdir: str
    values: dict[str, Value]


def save_job(job: Job, path: Path) -> None:
    payload: SavedJob = {"version": 1, "task": job.task, "workdir": str(job.workdir), "values": job.values}
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def load_job(path: Path) -> Job:
    path = path.expanduser().resolve()
    payload = cast(SavedJob, json.loads(path.read_text(encoding="utf-8")))
    if (not isinstance(payload, dict) or payload.get("version") != 1
            or not isinstance(payload.get("task"), str) or not isinstance(payload.get("workdir"), str)
            or not isinstance(payload.get("values"), dict)):
        raise ValueError("任务配置需要 version=1、task、workdir 和 values")
    workdir = Path(payload["workdir"]).expanduser()
    job = new_job(payload["task"], workdir if workdir.is_absolute() else path.parent / workdir)
    job.values.update(payload["values"])
    # Check schema without requiring that the user's data has been transferred yet.
    fields = {field.option: field for field in fields_for(job.task)}
    for option, value in job.values.items():
        if option not in fields or not isinstance(value, bool if fields[option].flag else str):
            raise ValueError(f"未知参数或值类型错误：{option}")
    return job


def run_job(job: Job) -> int:
    process = subprocess.Popen(build_command(job), cwd=job.workdir)
    try:
        return process.wait()
    except KeyboardInterrupt:
        if process.poll() is None:
            process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        return 130
