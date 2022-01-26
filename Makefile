.PHONY: clean All

All:
	@echo "----------Building project:[ static_lib_sv_rs_dopp - Debug ]----------"
	@cd "static_lib_sv_rs_dopp" && "$(MAKE)" -f  "static_lib_sv_rs_dopp.mk"
clean:
	@echo "----------Cleaning project:[ static_lib_sv_rs_dopp - Debug ]----------"
	@cd "static_lib_sv_rs_dopp" && "$(MAKE)" -f  "static_lib_sv_rs_dopp.mk" clean
