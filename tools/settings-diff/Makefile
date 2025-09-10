
PYUIC=pyuic5
PFLAGS=

all: ui/MainWindow_design.py

ui/%_design.py: ui/%.ui
	$(PYUIC) $(PFLAGS) $< -o $@

clean:
	rm ui/*.py
