"""Chinese terminal interface over the existing production entrypoints."""

from __future__ import annotations

import curses
from pathlib import Path
import unicodedata

from .tasks import (Job, STEP_TITLES, TASKS, Value, build_command, command_text,
                    fields_for, get_task, load_job, new_job, run_job, save_job)


def cell_width(text: str) -> int:
    return sum(0 if unicodedata.combining(char) else
               2 if unicodedata.east_asian_width(char) in {"W", "F"} else 1 for char in text)


def fit_text(text: str, cells: int) -> str:
    result: list[str] = []
    used = 0
    for char in text:
        char = char if char.isprintable() else " "
        width = cell_width(char)
        if used + width > cells:
            break
        result.append(char)
        used += width
    return "".join(result)


def wrap_text(text: str, cells: int) -> list[str]:
    lines: list[str] = []
    for paragraph in text.splitlines():
        if not paragraph:
            lines.append("")
        while paragraph:
            line = fit_text(paragraph, cells)
            if not line:
                break
            lines.append(line)
            paragraph = paragraph[len(line):]
    return lines


class QuitUI(Exception):
    pass


class TerminalApp:
    def __init__(self, screen: curses.window, workdir: Path, initial: Job | None = None) -> None:
        self.screen = screen
        self.workdir = workdir
        self.job = initial
        self.jobs: dict[str, Job] = {initial.task: initial} if initial else {}
        self.message = ""
        self.screen.keypad(True)
        curses.set_escdelay(30)
        self.cursor(False)
        self.accent = curses.A_BOLD
        if curses.has_colors():
            curses.start_color()
            curses.init_pair(1, curses.COLOR_CYAN, curses.COLOR_BLACK)
            self.accent |= curses.color_pair(1)

    def cursor(self, visible: bool) -> None:
        try:
            curses.curs_set(int(visible))
        except curses.error:
            pass

    def put(self, row: int, text: str, column: int = 2, style: int = 0) -> None:
        height, width = self.screen.getmaxyx()
        if not 0 <= row < height or not 0 <= column < width - 1:
            return
        try:
            self.screen.addstr(row, column, fit_text(text, width - column - 1), style)
        except curses.error:
            pass

    def ready(self) -> bool:
        while True:
            height, width = self.screen.getmaxyx()
            if height >= 16 and width >= 68:
                return True
            self.screen.erase()
            self.put(0, "请将终端扩大到至少 68 列 × 16 行；q 退出", 0)
            self.screen.refresh()
            if self.screen.get_wch() in {"q", "\x1b"}:
                raise QuitUI

    def frame(self, title: str) -> None:
        self.screen.erase()
        self.put(0, "CLE  |  统一生产入口", style=self.accent)
        self.put(2, title, style=curses.A_BOLD)

    def choose(self, title: str, labels: list[str]) -> int | None:
        index = 0
        while self.ready():
            self.frame(title)
            height, _ = self.screen.getmaxyx()
            rows = height - 7
            start = max(0, index - rows + 1)
            for row, label in enumerate(labels[start:start + rows], start=4):
                selected = row - 4 + start == index
                self.put(row, f"{'›' if selected else ' '} {label}", style=curses.A_REVERSE if selected else 0)
            self.put(height - 2, "↑↓ / j k 选择    Enter 进入    Esc 返回    q 退出")
            self.screen.refresh()
            key = self.screen.get_wch()
            if key in {curses.KEY_UP, "k"}:
                index = max(0, index - 1)
            elif key in {curses.KEY_DOWN, "j"}:
                index = min(len(labels) - 1, index + 1)
            elif key in {"\n", "\r", curses.KEY_ENTER}:
                return index
            elif key == "\x1b":
                return None
            elif key == "q":
                raise QuitUI
        return None

    def edit(self, title: str, initial: str, hint: str = "") -> str | None:
        text, position, start = initial, len(initial), 0
        self.cursor(True)
        try:
            while self.ready():
                self.frame(title)
                height, width = self.screen.getmaxyx()
                self.put(4, hint)
                available = width - 5
                start = min(start, position)
                while cell_width(text[start:position]) >= available:
                    start += 1
                self.put(7, fit_text(text[start:], available), style=curses.A_UNDERLINE)
                self.put(height - 3, "Enter 保存    Esc 取消    Ctrl+U 清空    Home/End 或 Ctrl+A/E")
                self.screen.move(7, 2 + cell_width(text[start:position]))
                self.screen.refresh()
                key = self.screen.get_wch()
                if key in {"\n", "\r", curses.KEY_ENTER}:
                    return text
                if key == "\x1b":
                    return None
                if key in {curses.KEY_BACKSPACE, "\x7f", "\b"}:
                    if position:
                        text = text[:position - 1] + text[position:]
                        position -= 1
                elif key == curses.KEY_DC:
                    text = text[:position] + text[position + 1:]
                elif key == curses.KEY_LEFT:
                    position = max(0, position - 1)
                elif key == curses.KEY_RIGHT:
                    position = min(len(text), position + 1)
                elif key in {curses.KEY_HOME, "\x01"}:
                    position = 0
                elif key in {curses.KEY_END, "\x05"}:
                    position = len(text)
                elif key == "\x15":
                    text, position, start = "", 0, 0
                elif isinstance(key, str) and key.isprintable():
                    text = text[:position] + key + text[position:]
                    position += 1
        finally:
            self.cursor(False)
        return None

    def config_path(self, text: str) -> Path:
        path = Path(text).expanduser()
        workdir = self.job.workdir if self.job else self.workdir
        return path if path.is_absolute() else workdir / path

    def save(self) -> None:
        assert self.job is not None
        default = self.job.workdir / "work" / "tui" / f"{self.job.task}.json"
        text = self.edit("保存任务配置（使用新文件名）", str(default), "保存当前参数与工作目录，不启动计算")
        if text is not None and text.strip():
            try:
                path = self.config_path(text)
                save_job(self.job, path)
                self.message = f"已保存：{path}"
            except (OSError, ValueError) as error:
                self.message = f"保存失败：{error}"

    def load(self) -> None:
        text = self.edit("载入任务配置", "", "输入任务配置 JSON 路径；相对路径以工作目录为准")
        if text is not None and text.strip():
            try:
                job = load_job(self.config_path(text))
                self.jobs[job.task] = job
                self.job = job
                self.message = f"已载入：{job.task}"
            except (OSError, ValueError) as error:
                self.message = f"载入失败：{error}"

    def execute(self) -> None:
        assert self.job is not None
        command = build_command(self.job)
        curses.def_prog_mode()
        curses.endwin()
        try:
            print(command_text(self.job, command), flush=True)
            code = run_job(self.job)
            self.message = "任务已中断（退出码 130）" if code == 130 else f"任务已结束，退出码 {code}"
            print(f"\n{self.message}", flush=True)
            try:
                input("按 Enter 返回参数表单：")
            except (KeyboardInterrupt, EOFError):
                pass
        except (OSError, ValueError) as error:
            self.message = f"启动失败：{error}"
        finally:
            curses.reset_prog_mode()
            self.screen.keypad(True)
            self.screen.clear()

    def preview(self) -> None:
        assert self.job is not None
        try:
            command = build_command(self.job)
        except ValueError as error:
            self.message = str(error)
            return
        offset = 0
        while self.ready():
            self.frame(f"命令预览  |  {get_task(self.job.task).title}")
            height, width = self.screen.getmaxyx()
            lines = wrap_text(command_text(self.job, command), width - 5)
            rows = height - 8
            offset = min(offset, max(0, len(lines) - rows))
            for row, line in enumerate(lines[offset:offset + rows], start=4):
                self.put(row, line)
            self.put(height - 3, "r 启动此命令    s 保存参数    Enter / Esc 返回修改")
            self.put(height - 2, "↑↓ / j k 滚动    q 退出；运行后原脚本输出直接显示在终端")
            self.screen.refresh()
            key = self.screen.get_wch()
            if key in {"\n", "\r", "\x1b", curses.KEY_ENTER}:
                return
            if key == "r":
                try:
                    self.execute()
                except (OSError, ValueError) as error:
                    self.message = f"启动失败：{error}"
                return
            if key == "s":
                self.save()
            elif key in {curses.KEY_DOWN, "j"}:
                offset = min(max(0, len(lines) - rows), offset + 1)
            elif key in {curses.KEY_UP, "k"}:
                offset = max(0, offset - 1)
            elif key == "q":
                raise QuitUI

    def form(self) -> None:
        index, start = 0, 0
        while self.ready():
            assert self.job is not None
            fields = fields_for(self.job.task)
            self.frame(f"{get_task(self.job.task).title}  |  {self.job.task}.py")
            height, width = self.screen.getmaxyx()
            entries = [("工作目录", str(self.job.workdir))]
            for field in fields:
                value: Value = self.job.values.get(field.option, field.initial)
                displayed = ("[x]" if value else "[ ]") if field.flag else str(value) or "〈未填写〉"
                entries.append((f"{'*' if field.action.required else ' '} {field.label}", displayed))
            index = min(index, len(entries) - 1)
            rows = height - 10
            start = max(min(start, index), index - rows + 1)
            label_width = min(34, width // 2)
            for row, (label, value) in enumerate(entries[start:start + rows], start=4):
                selected = row - 4 + start == index
                style = curses.A_REVERSE if selected else 0
                label = fit_text(label, label_width - 2)
                self.put(row, label + " " * (label_width - cell_width(label)) + value, style=style)
            hint = fields[index - 1].hint if index else "工作目录决定相对数据路径的解释方式"
            self.put(height - 6, hint)
            for row, line in enumerate(wrap_text(self.message, width - 5)[:2], start=height - 5):
                self.put(row, line)
            self.put(height - 2, "↑↓ / j k 选择  Enter 编辑  p / r 命令预览  s 保存  l 载入")
            self.put(height - 1, "Esc 返回步骤选择  q 退出  * 必填；空白可选项使用脚本默认值")
            self.screen.refresh()
            key = self.screen.get_wch()
            if key in {curses.KEY_UP, "k"}:
                index = max(0, index - 1)
            elif key in {curses.KEY_DOWN, "j", "\t"}:
                index = min(len(entries) - 1, index + 1)
            elif key == curses.KEY_NPAGE:
                index = min(len(entries) - 1, index + rows)
            elif key == curses.KEY_PPAGE:
                index = max(0, index - rows)
            elif key in {"\n", "\r", curses.KEY_ENTER, " "}:
                self.message = ""
                if index == 0:
                    text = self.edit("工作目录", str(self.job.workdir))
                    if text is not None and text.strip():
                        self.job.workdir = Path(text).expanduser().resolve()
                else:
                    field = fields[index - 1]
                    value = self.job.values.get(field.option, field.initial)
                    if field.flag:
                        self.job.values[field.option] = not value
                    else:
                        text = self.edit(f"{field.label}  |  {field.option}", str(value), field.hint)
                        if text is not None:
                            self.job.values[field.option] = text
            elif key in {"p", "r"}:
                self.preview()
            elif key == "s":
                self.save()
            elif key == "l":
                self.load()
                index, start = 0, 0
            elif key == "\x1b":
                self.job = None
                return
            elif key == "q":
                raise QuitUI

    def run(self) -> None:
        try:
            while True:
                if self.job is None:
                    steps = list(STEP_TITLES)
                    selected = self.choose("选择生产步骤", list(STEP_TITLES.values()))
                    if selected is None:
                        return
                    tasks = [task for task in TASKS if task.step == steps[selected]]
                    chosen = self.choose(STEP_TITLES[steps[selected]], [f"{task.title}  ({task.name})" for task in tasks])
                    if chosen is None:
                        continue
                    name = tasks[chosen].name
                    if name not in self.jobs:
                        self.jobs[name] = new_job(name, self.workdir)
                    self.job = self.jobs[name]
                    self.message = ""
                self.form()
        except QuitUI:
            return
