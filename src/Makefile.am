SUBDIRS = SGTpsv SGTsh KernelMaker ModelDrawer ReadKernel
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm $(top_builddir)/bin/*.optrpt; \ 
	for dir in $(SUBDIRS); do \
	rm $$dir/*.optrpt; \
	rm $$dir/*.mod; \
	done  
