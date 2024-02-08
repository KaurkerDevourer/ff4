SUBDIRS = 	buchberger \
			rational

all: $(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

run: 
	./bin/test_buchberger && \
	./bin/test_rational

.PHONY: all $(SUBDIRS) run