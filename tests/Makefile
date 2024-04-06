SUBDIRS = 	buchberger \
			rational \
			term \
			monomial \
			polynomial

all: $(SUBDIRS)

$(SUBDIRS):
		$(MAKE) -C $@

run: 
	./bin/test_rational && \
	./bin/test_term && \
	./bin/test_monomial && \
	./bin/test_polynomial && \
	./bin/test_buchberger

.PHONY: all $(SUBDIRS) run
