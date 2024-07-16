import pytest

from ElliptiCBn.run_all import _file_check
from ElliptiCBn.run_all import run_all

import os

def test__file_check(tmpdir):
    
    cwd = os.getcwd()
    os.chdir(tmpdir)

    output = _file_check(some_file="no_path_file.txt",
                         output_dir="output_dir")
    assert output == os.path.join("output_dir","no_path_file.txt")

    output = _file_check(some_file="../with/path/file.txt",
                         output_dir="output_dir")
    assert output == "../with/path/file.txt"

    output = _file_check(some_file="../with/path/file.txt",
                         output_dir="something_else")
    assert output == "../with/path/file.txt"

    with open("file-exists.txt","w") as f:
        f.write("line\n")
    assert os.path.isfile("file-exists.txt")

    with pytest.raises(FileExistsError):
        output = _file_check(some_file="./file-exists.txt")

    with pytest.raises(FileExistsError):
        output = _file_check(some_file="file-exists.txt",
                             output_dir=".")
        
    output = _file_check(some_file="./file-exists.txt",overwrite=True)
    assert output == "./file-exists.txt"
    assert not os.path.isfile("file-exists.txt")

    os.chdir(cwd)


def test_run_all(example_xyz,tmpdir):
    
    cwd = os.getcwd()
    os.chdir(tmpdir)

    for counter, xyz in enumerate(example_xyz["*.xyz"]):

        # no file exists already
        out_file = os.path.join(".",f"{os.path.split(xyz)[-1]}.html")
        assert not os.path.isfile(out_file)

        # should work
        run_all(xyz,
                output_dir=".",
                overwrite=False)
        assert os.path.isfile(out_file)

        # should fail -- we just made output
        with pytest.raises(FileExistsError):
            run_all(xyz,
                    output_dir=".",
                    overwrite=False)
            
        # should work now -- overwrite
        run_all(xyz,
                output_dir=".",
                overwrite=True)

        with pytest.raises(ValueError):
            run_all(xyz,
                    min_num_carbons=10000,
                    max_num_carbons=10,
                    output_dir=".",
                    overwrite=True)
            
    # Make sure it writes out blocks of files to directory with summary
    assert not os.path.isdir("output-files")
    run_all(example_xyz["*.xyz"],output_dir="output-files",summary_file="custom.xlsx")
    assert os.path.isdir("output-files")
    assert os.path.isfile(os.path.join("output-files","custom.xlsx"))

    # summary csv
    run_all(example_xyz["*.xyz"],
            output_dir="output-files",
            summary_file="custom.csv",
            overwrite=True)
    assert os.path.isfile(os.path.join("output-files","custom.csv"))

    # Fail -- can't figure out csv or xlsx
    with pytest.raises(ValueError):
        run_all(example_xyz["*.xyz"],
                output_dir="output-files",
                summary_file="custom.blah",
                overwrite=True)
        

        
    with open("not-a-dir.txt","w") as f:
        f.write('here')

    with pytest.raises(FileExistsError):
        run_all(example_xyz["*.xyz"],
                output_dir="not-a-dir.txt",
                overwrite=False)
        
    # Wipe out file and replace with directory
    assert os.path.isfile("not-a-dir.txt")
    assert not os.path.isdir("not-a-dir.txt")
    run_all(example_xyz["*.xyz"],
            output_dir="not-a-dir.txt",
            overwrite=True)
    assert os.path.isdir("not-a-dir.txt")

    with pytest.raises(ValueError):
        run_all(None,
                output_dir="not-a-dir.txt",
                overwrite=True)


    os.chdir(cwd)
