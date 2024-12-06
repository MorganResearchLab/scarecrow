import logging
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
    logger.setLevel(log_level)

    # Clear any existing handlers to prevent duplicate logging
    logger.handlers.clear()

    # File Handler
    try:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(log_level)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        
        # Add file handler to logger
        logger.addHandler(file_handler)
        
        # Add a stream handler to capture any potential logging errors
        stream_handler = logging.StreamHandler(sys.stderr)
        stream_handler.setLevel(logging.ERROR)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
        
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
logger = setup_logger()

# Demonstration of logging usage
@log_errors
def example_function():
    logger.info("This is an informational message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    
    # Simulate an error
    x = 1 / 0

if __name__ == "__main__":
    try:
        example_function()
    except Exception:
        print("Function execution failed. Check log file.")