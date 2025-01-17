import logging
from rich.logging import RichHandler
from logging.handlers import WatchedFileHandler  # Add this import
import functools
import os
import sys

def setup_logger(log_file: str = 'scarecrow.log',
    log_level: int = logging.INFO
) -> logging.Logger:
    """
    Create a logger that writes to file with robust configuration.
    
    Args:
        log_file (str): Path to log file
        log_level (int): Logging level
    
    Returns:
        logging.Logger: Configured logger
    """
    # Ensure log directory exists
    log_dir = os.path.dirname(log_file) or '.'
    os.makedirs(log_dir, exist_ok=True)

    # Create logger
    logger = logging.getLogger('scarecrow')

    # If logger already has handlers, remove them
    if logger.hasHandlers():
        logger.handlers.clear()

    try:
        shell_handler = RichHandler()
        # Use WatchedFileHandler instead of FileHandler for better multiprocess support
        file_handler = WatchedFileHandler(log_file, mode='a')
        logger.setLevel(log_level)
        shell_handler.setLevel(log_level)
        file_handler.setLevel(log_level)

        # Clear any existing handlers to prevent duplicate logging
        logger.handlers.clear()
        
        # Create formatter
        fmt_shell = '%(message)s'
        fmt_file = '%(levelname)s %(asctime)s [%(filename)s:%(funcName)s:%(lineno)d] %(message)s'
        shell_formatter = logging.Formatter(fmt_shell)
        file_formatter = logging.Formatter(fmt_file)
        shell_handler.setFormatter(shell_formatter)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(shell_handler)
        logger.addHandler(file_handler)

        # Add log file name to logger
        logger.filename = log_file
 
        # Test logging
        logger.info("Logger successfully initialized")

    except Exception as e:
        # Fallback error handling
        print(f"CRITICAL: Unable to set up logger: {e}", file=sys.stderr)
        raise

    return logger

# Example usage and error handling decorator
def log_errors(func):
    """
    Decorator to log any errors that occur in the decorated function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Error in {func.__name__}: {e}", exc_info=True)
            raise
    return wrapper


# Global
#logger = setup_logger()
