from pathlib import Path
import subprocess
import sys


script = Path("01_pipeline/scripts/11_v14_qc_audit.py")
cmd = [sys.executable, str(script), *sys.argv[1:]]
raise SystemExit(subprocess.run(cmd).returncode)
