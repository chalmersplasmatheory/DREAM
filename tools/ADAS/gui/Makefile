
PYUIC=pyuic5
PFLAGS=

all: ManualFit_design.py

%_design.py: %.ui
	$(PYUIC) $(PFLAGS) $< -o $@

clean:
	rm *_design.py

