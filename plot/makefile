CFLAGS = -w
LDFLAGS = -lm -lpthread
generate_cores_1:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o one_core ${argument} $(LDFLAGS)
generate_cores_2:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o two_core ${argument} $(LDFLAGS)
generate_cores_3:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o third_core ${argument} $(LDFLAGS)
generate_cores_4:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o fourth_core ${argument} $(LDFLAGS)
generate_cores_5:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o fifth_core ${argument} $(LDFLAGS)
generate_cores_6:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o six_core ${argument} $(LDFLAGS)
generate_cores_7:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o seven_core ${argument} $(LDFLAGS)
generate_cores_8:${argument}
	@gcc $(CFLAGS) -O3 -ffast-math -msse4.2 -o eigth_core ${argument} $(LDFLAGS)


clean:
	@rm -f one_core
	@rm -f two_core
	@rm -f third_core
	@rm -f fourth_core
	@rm -f fifth_core
	@rm -f six_core
	@rm -f seven_core
	@rm -f eigth_core
	@rm time.txt
