
PYUIC=pyuic5
PFLAGS=

all: ui/EqView_design.py ui/EqViewConfig_design.py ui/MainWindow_design.py

ui/%_design.py: ui/%.ui
	$(PYUIC) $(PFLAGS) $< -o $@

clean:
	rm ui/*.py

