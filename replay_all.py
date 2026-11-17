"""Replay pinned finite certificates; analytic hypotheses remain external."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent
    output = Path(os.environ.get('SM_OUT', str(root / 'logs' / 'merged-replay'))).resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for name, expected in json.loads((root / 'INPUT-HASHES.json').read_text()).items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit(f'INPUT HASH MISMATCH: {name}')
    print('PINNED INPUT HASHES PASS', flush=True)
    degree4 = root / 'scripts'
    if not (degree4 / 'cert_barriers.py').is_file():
        degree4 = root / 'certificates'  # portable release layout
    commands = [(degree4, [name, '4']) for name in
                ['cert_barriers.py', 'cert_upper.py', 'cert_lower.py', 'sanity_checks.py']]
    commands += [(degree4, ['verify_flint.py']),
                 (root / 'degrees5-8', ['verify_flint.py', 'certificates']),
                 (root / 'degrees5-8', ['semantic_mutants.py'])]
    for cwd, argv in commands:
        cmd = [sys.executable, '-B', *argv]
        print('RUN', cwd.relative_to(root), *argv, flush=True)
        env = dict(os.environ, SM_OUT=str(output))
        result = subprocess.run(cmd, cwd=cwd, env=env)
        records.append({'cwd': cwd.relative_to(root).as_posix(), 'argv': argv, 'exit': result.returncode})
        (output / 'exit-codes.json').write_text(json.dumps(records, indent=2) + '\n')
        if result.returncode:
            return result.returncode
    print('EXACT FINITE REPLAY PASS; analytic hypotheses and variational bridge remain external.', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
