
import pytest

import os
import subprocess

def test_script(example_xyz,tmpdir):
    """
    Basic test of the command line script.
    """

    cwd = os.getcwd()
    os.chdir(tmpdir)

    expected = []
    for xyz in example_xyz["*.xyz"]:

        out_html = f'{os.path.split(xyz)[-1]}.html'
        out_xlsx = f'{os.path.split(xyz)[-1]}.xlsx'

        expected.append(out_html)
        expected.append(out_xlsx)

        cmd = ["ElliptiCBn",xyz]
        results = subprocess.run(cmd)
        assert results.returncode == 0
        assert os.path.isfile(out_html)
        assert os.path.isfile(out_xlsx)

        os.remove(out_html)
        os.remove(out_xlsx)

    cmd = ["ElliptiCBn",*example_xyz["*.xyz"]]
    results = subprocess.run(cmd)
    assert results.returncode == 0

    expected.append("summary.xlsx")
    assert set(os.listdir(".")) == set(expected)
    
    # make sure we can write to a directory
    cmd = ["ElliptiCBn",*example_xyz["*.xyz"]]
    cmd.extend(["--output_dir","multi"])
    results = subprocess.run(cmd)
    assert results.returncode == 0
    assert os.path.isdir("multi")
    assert set(os.listdir("multi")) == set(expected)

    # Should fail -- output exists
    cmd = ["ElliptiCBn",*example_xyz["*.xyz"]]
    cmd.extend(["--output_dir","multi"])
    results = subprocess.run(cmd)
    assert results.returncode == 1

    # Should work now
    cmd.append("--overwrite")
    results = subprocess.run(cmd)
    assert results.returncode == 0
    

    cmd = ["ElliptiCBn",example_xyz["*.xyz"][0]]
    cmd.extend(["--min_num_carbons","1000"])
    cmd.extend(["--max_num_carbons","20"])
    results = subprocess.run(cmd)
    assert results.returncode == 1

    cmd = ["ElliptiCBn",example_xyz["*.xyz"][0]]
    cmd.extend(["--min_num_carbons","20"])
    cmd.extend(["--max_num_carbons","2"])
    results = subprocess.run(cmd)
    assert results.returncode == 1




    os.chdir(cwd)


