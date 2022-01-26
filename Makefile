.PHONY: clean All

All:
	@echo "----------Building project:[ lib_sv_rs_dopp - Debug ]----------"
	@cd "lib_sv_rs_dopp" && "$(MAKE)" -f  "lib_sv_rs_dopp.mk"
clean:
	@echo "----------Cleaning project:[ lib_sv_rs_dopp - Debug ]----------"
	@cd "lib_sv_rs_dopp" && "$(MAKE)" -f  "lib_sv_rs_dopp.mk" clean
