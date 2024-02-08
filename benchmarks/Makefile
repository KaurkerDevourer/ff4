SUBDIRS = 	buchberger \
			sum

all: $(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

run: 
	./bin/benchmark_buchberger && \
	./bin/benchmark_sum

.PHONY: all $(SUBDIRS) run