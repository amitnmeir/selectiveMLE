.PHONY: install clean

install: clean
	@R CMD INSTALL .

clean:
	@echo "Cleaning build artifacts"
	rm -f src/*.o src/*.so src/*.dll
	@echo "Done."
