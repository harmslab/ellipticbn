
import pytest

import os
import subprocess
import platform


def test_script(example_xyz,tmpdir):
    """
    Basic test of the command line script.
    """

    # # Update test for windows
    # if platform.system() == "Windows":
    #     here = os.path.abspath(os.path.dirname(__file__))
    #     script = os.path.realpath(os.path.join(here,"..","bin","ElliptiCBn"))
    #     base_cmd = ["python",script]
    # else:
    base_cmd = ["ElliptiCBn"]

    cwd = os.getcwd()
    os.chdir(tmpdir)

    expected = []
    for xyz in example_xyz["*.xyz"]:

        out_html = f'{os.path.split(xyz)[-1]}.html'
        out_xlsx = f'{os.path.split(xyz)[-1]}.xlsx'

        expected.append(out_html)
        expected.append(out_xlsx)

        cmd = base_cmd[:]
        cmd.append(xyz)
        results = subprocess.run(cmd)
        assert results.returncode == 0
        assert os.path.isfile(out_html)
        assert os.path.isfile(out_xlsx)

        os.remove(out_html)
        os.remove(out_xlsx)

    cmd = base_cmd[:]
    cmd.extend(example_xyz["*.xyz"])
    results = subprocess.run(cmd)
    assert results.returncode == 0

    expected.append("summary.xlsx")
    seen = set(os.listdir("."))
    seen = seen - set(["tmp-xyz.xyz"])
    assert seen == set(expected)
    
    # make sure we can write to a directory
    cmd = base_cmd[:]
    cmd.extend(example_xyz["*.xyz"])
    cmd.extend(["--output_dir","multi"])
    results = subprocess.run(cmd)
    assert results.returncode == 0
    assert os.path.isdir("multi")

    # Make sure the right files are spit out (less a temporary file that might
    # not be deleted)
    seen = set(os.listdir("multi"))
    seen = seen - set(["tmp-xyz.xyz"])
    assert seen == set(expected)

    # Should fail -- output exists
    cmd = base_cmd[:]
    cmd.extend(example_xyz["*.xyz"])
    cmd.extend(["--output_dir","multi"])
    results = subprocess.run(cmd)
    assert results.returncode == 1

    # Should work now
    cmd.append("--overwrite")
    results = subprocess.run(cmd)
    assert results.returncode == 0
    
    cmd = base_cmd[:]
    cmd.append(example_xyz["*.xyz"][0])
    cmd.extend(["--min_num_carbons","1000"])
    cmd.extend(["--max_num_carbons","20"])
    results = subprocess.run(cmd)
    assert results.returncode == 1

    cmd = base_cmd[:]
    cmd.append(example_xyz["*.xyz"][0])
    cmd.extend(["--min_num_carbons","20"])
    cmd.extend(["--max_num_carbons","2"])
    results = subprocess.run(cmd)
    assert results.returncode == 1

    os.chdir(cwd)


