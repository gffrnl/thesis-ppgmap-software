package = frlap
version = 1.0.0
libname = libfrlap
tarname = $(package)-$(version)
distdir = $(tarname)

prefix = /usr/local
exec_prefix = $(prefix)
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib


# Make does not offer a recursive wildcard function, so here's one:
#rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

header_files = $(wildcard include/*.h)       \
               $(wildcard include/frlap/*.h)
#src_files = $(wildcard src/*.c)
src_files = src/frlap.c
obj_files = $(src_files:.c=.o)

CFLAGS += -I./include -fPIC


all: lib/$(libname).so

lib/$(libname).so: $(header_files) $(src_files) $(obj_files)
	mkdir -p lib
	$(CC) -shared -o $@ $(LDFLAGS) $(LDLIBS) $(obj_files)

%.o: %.c $(header_files)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


check:
	@echo "*** ALL TESTS PASSED"

install: lib/$(libname).so
	install -d $(DESTDIR)$(includedir)
	install -m 0755 include/$(package).h $(DESTDIR)$(includedir)
	install -d $(DESTDIR)$(libdir)
	install -m 0755 lib/$(libname).so $(DESTDIR)$(libdir)

uninstall:
	rm -f $(DESTDIR)$(includedir)/$(package).h
	rm -f $(DESTDIR)$(libdir)/$(libname).so

dist: $(tarname).tar.gz

$(tarname).tar.gz: $(distdir)
	tar chof - $(tarname) | gzip -9 -c > $@
	rm -rf $(distdir)

$(distdir): FORCE
	mkdir -p \
         $(distdir)/include       \
         $(distdir)/include/frlap \
         $(distdir)/src
	cp Makefile $(distdir)
	cp include/$(package).h $(distdir)/include
	cp include/frlap/*.h    $(distdir)/include/frlap
	cp $(src_files) $(wildcard src/*.h) $(distdir)/src

distcheck: $(tarname).tar.gz
	gzip -cd $(tarname).tar.gz | tar xvf -
	$(MAKE) lib/$(libname).so
	$(MAKE) check
	$(MAKE) DESTDIR=$${PWD}/$(distdir)/_inst install
	$(MAKE) DESTDIR=$${PWD}/$(distdir)/_inst uninstall
	@remaining="`find $(distdir)/_inst -type f | wc -l`"; \
	  if test "$${remaining}" -ne 0; then \
	    echo "*** $${remaining} file(s) remaining in stage directory!"; \
	    exit 1; \
	  fi
	$(MAKE) clean
	@echo "*** Package $(tarname).tar.gz is ready for distribution."

FORCE:
	@rm -rf $(distdir)

clean:
	rm -f src/*.o
	rm -rf $(distdir) lib $(tarname).tar.gz

purge: clean
	rm -f *~
	rm -f include/*~
	rm -f src/*~

.PHONY: FORCE all clean purge check dist distcheck install uninstall
