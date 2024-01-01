# Define directories
DIRS := cdt cgraph circogen common dotgen fdpgen gvc label \
neatogen ortho osage pack patchwork pathplan rbtree sfdpgen sparse \
twopigen vpsc xdot \
plugin_core plugin_dot_layout plugin_neato_layout

all: $(DIRS)

# Pattern rule for directories
$(DIRS):
	$(MAKE) -C $@ all

.PHONY: all $(DIRS)

# Clean rule
clean:
	for dir in $(DIRS); do \
		$(MAKE) -C $$dir clean; \
	done