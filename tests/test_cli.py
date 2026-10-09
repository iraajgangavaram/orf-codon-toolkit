from orfkit.cli import main


def _fasta(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">rec1 test\nGGATGAAAAAAAAATAA\n")
    return str(p)


def test_orfs_command(tmp_path, capsys):
    assert main(["orfs", _fasta(tmp_path), "--min-aa", "4", "--forward-only"]) == 0
    rows = capsys.readouterr().out.strip().split("\n")
    assert rows[0].startswith("record\tstrand")
    assert rows[1].split("\t") == ["rec1", "+", "3", "3", "17", "4", "MKKK"]


def test_gc_command(tmp_path, capsys):
    main(["gc", _fasta(tmp_path)])
    assert capsys.readouterr().out.strip().split("\n")[1] == "rec1\t17\t0.1765"


def test_broken_pipe_is_silent():
    import os
    import subprocess
    import sys

    root = os.path.dirname(os.path.dirname(__file__))
    fasta = os.path.join(root, "examples", "example.fa")
    cmd = f"{sys.executable} -m orfkit codons {fasta} | head -n 1"
    r = subprocess.run(cmd, shell=True, cwd=root, capture_output=True, text=True)
    assert "Traceback" not in r.stderr and r.stdout.startswith("record")
