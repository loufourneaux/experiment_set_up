import os
import logging

def log_parser_isfile(file_path, parser, isForce):
    """
    Check if the file exists at the given file path and handle accordingly.

    Parameters:
    file_path (str): The path to the file.
    parser (argparse.ArgumentParser): The argument parser instance.
    isForce (bool): Flag to force overwrite of the file if it exists.

    Raises:
    SystemExit: If the file exists and isForce is False, raising a parser error.
    """
    if os.path.isfile(file_path):
        if isForce:
            logging.info('Overwriting "{0}".'.format(file_path))
        else:
            parser.error(
                '"{0}" already exists! Use -f to overwrite it.'
                .format(file_path))
