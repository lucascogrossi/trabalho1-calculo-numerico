main: main.c tinyexpr/tinyexpr.c
	@gcc main.c tinyexpr/tinyexpr.c -o main -lm

clean:
	@rm -f main

run: main
	@./main

