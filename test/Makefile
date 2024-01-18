
help:
	@echo tests

# run installed dot in order to generate comparison
setup:
	dot -Tjpg -o arrows-good.jpg arrows.gv
	dot -Tjpg -o b-good.jpg b.gv
	open *.jpg

test:
	./devdot.sh -Tsvg -o arrows.svg arrows.gv
	./devdot.sh -Tsvg -o b.svg b.gv
	open *.svg

clean:
	*.svg
