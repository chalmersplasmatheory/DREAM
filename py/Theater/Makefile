
PYUIC=pyuic5
PFLAGS=

all: ui/DREAMTheater_design.py

ui/%_design.py: ui/%.ui
	$(PYUIC) $(PFLAGS) $< -o $@

clean:
	rm ui/*_design.py

