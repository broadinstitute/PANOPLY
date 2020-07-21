## install   : install tedmint
.PHONY: install
install:
	bash install_tool.sh

## uninstall : uninstall tedmint
.PHONY: uninstall
uninstall:
	rm -rf ~/.tedmint

## help      : display help message
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<
