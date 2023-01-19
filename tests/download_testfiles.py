import sys
import directlfq.testfile_handling
test_folder = "."
links_yaml_all_testfiles = "../directlfq/configs/download_links_for_testfiles_all.yaml"
links_yaml_quicktest_files = "../directlfq/configs/download_links_for_testfiles_quicktest.yaml"

def get_links_yaml_file_as_specified(command_line_arguments):
    if command_line_arguments[1] == 'quicktest':
        return links_yaml_quicktest_files
    if command_line_arguments[1] == 'all_tests':
        return links_yaml_all_testfiles


if __name__ == '__main__':
    command_line_arguments = sys.argv
    try:
        links_yaml = get_links_yaml_file_as_specified(command_line_arguments)
    except:
        raise ValueError("specify if \'quicktest\' or \'all_tests\' on command line")

    testfieldownloader = directlfq.testfile_handling.TestFolderDownloader(test_folder=test_folder, links_yaml=links_yaml)
    testfieldownloader.download_missing_files()