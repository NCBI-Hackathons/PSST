all: build

build:
	g++ -std=c++11 -o unit_tests catch_config_main.cpp test_data_manipulation.cpp test_queries_with_ref_bases.cpp test_variant_calling.cpp ../data_manipulation.cpp ../queries_with_ref_bases.cpp ../variant_calling.cpp

clean:
	rm *.o unit_tests
