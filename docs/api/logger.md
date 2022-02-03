# The `logger` submodule

## Custom logging classes

::: safedata_validator.logger
    rendering:
        show_root_heading: false
        show_source: true
    selection:
        members:
            - CounterHandler
            - IndentFormatter

## Global logging instances

::: safedata_validator.logger
    rendering:
        show_source: true
    selection:
        members:
            - CH
            - CONSOLE_LOG
            - LOG
            - FORMATTER

## Convenience functions

::: safedata_validator.logger
    rendering:
        show_source: true
    selection:
        members:
            - log_and_raise
            - loggerinfo_push_pop
