"""环肽、受体准备与对接的统一终端界面。"""

from __future__ import annotations

import argparse
import locale
import os
from pathlib import Path
import sys

from production.tasks import STEP_TITLES, TASKS, build_command, command_text, load_job, new_job


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=[task.name for task in TASKS], help="直接打开指定任务表单")
    parser.add_argument("--workdir", type=Path, help="数据路径的工作目录，默认当前目录或配置中的目录")
    parser.add_argument("--config", type=Path, help="载入已保存的任务配置 JSON")
    parser.add_argument("--list", action="store_true", help="列出可用任务，不打开界面")
    parser.add_argument("--print-command", action="store_true", help="打印配置对应的命令，不执行")
    args = parser.parse_args()
    if args.list:
        for step, title in STEP_TITLES.items():
            print(title)
            for task in TASKS:
                if task.step == step:
                    print(f"  {task.name:<24} {task.title}")
        return
    try:
        job = load_job(args.config) if args.config else None
        if job and args.task and job.task != args.task:
            raise ValueError("--task 与配置中的任务不一致")
        workdir = (args.workdir or (job.workdir if job else Path.cwd())).expanduser().resolve()
        if job:
            job.workdir = workdir
        elif args.task:
            job = new_job(args.task, workdir)
        if args.print_command:
            if job is None:
                raise ValueError("--print-command 需要 --config 或 --task")
            print(command_text(job, build_command(job, check_paths=False)))
            return
        if not sys.stdin.isatty() or not sys.stdout.isatty():
            raise ValueError("TUI 需要交互终端；Conda 请使用 conda run --no-capture-output -n cle python src/tui.py")
        if os.environ.get("TERM", "dumb") in {"dumb", "unknown"}:
            raise ValueError("终端没有可用的 TERM 设置，请在普通终端或 SSH 终端启动")
        import curses
        from production.tui import TerminalApp

        locale.setlocale(locale.LC_ALL, "")
        curses.wrapper(lambda screen: TerminalApp(screen, workdir, job).run())
    except (OSError, ValueError, ImportError) as error:
        parser.exit(1, f"错误：{error}\n")
    except KeyboardInterrupt:
        print("\n已退出终端界面。")


if __name__ == "__main__":
    main()
