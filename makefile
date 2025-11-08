projects = phase_common phase_rare switch ligate simulate xcftools tests

.PHONY: all $(projects) coverage

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done
	rm -f tests/coverage/*
	rm -f static_bins/*
	rm -f docker/resources/*
	rm -f docker/shapeit5*.tar.gz

coverage: $(addsuffix .cov,$(projects))

# each of these is an independent job; -j 24 will run many at once
%.cov:
	$(MAKE) -C $* COVERAGE=1


static_exe:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

