"""Run an explicit command with separate immutable logs for each attempt."""
import datetime
import json
import pathlib
import subprocess
import sys
import time

directory = pathlib.Path(sys.argv[1])
directory.mkdir(parents=True, exist_ok=False)
command = sys.argv[2:]
record = {"command": command, "cwd": str(pathlib.Path.cwd()),
          "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
(directory / "status.json").write_text(json.dumps(record, indent=2))
with (directory / "stdout.log").open("w") as out, (directory / "stderr.log").open("w") as err:
    process = subprocess.Popen(command, stdout=out, stderr=err)
    record["pid"] = process.pid
    (directory / "status.json").write_text(json.dumps(record, indent=2))
    result = process.wait()
record.update(exit_code=result, ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(directory / "status.json").write_text(json.dumps(record, indent=2))
print(json.dumps(record, indent=2))
print((directory / "stdout.log").read_text()[-6000:])
print((directory / "stderr.log").read_text()[-6000:])
sys.exit(result)
