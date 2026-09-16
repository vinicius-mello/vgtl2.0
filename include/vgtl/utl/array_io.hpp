#ifndef VGTL_ARRAY_IO_HPP
#define VGTL_ARRAY_IO_HPP

/*! \file
 * \brief Array IO operations
 */

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <vgtl/comb/range.hpp>
#include <fstream>

namespace vgtl {

	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::size_t;
	using boost::array;
	using boost::multi_array;
	
	/*! \addtogroup utl 
	 * @{
	 */

	//! Reads an array from a stream
	template <size_t k, class Data>
	void read(ifstream& in, array<Data,k>& a) {
		for(int i=0;i<k;++i) {
			in.read(reinterpret_cast<char *>(&a[i]),sizeof(Data));
		} 
	}
	
	//! Writes an array to a stream
	template <size_t k, class Data>
	void write(ofstream& out, array<Data,k>& a) {
		for(int i=0;i<k;++i) {
			out.write(reinterpret_cast<char *>(&a[i]),sizeof(Data));
		} 
	}
	
	//! Reads a multi array from a stream
	template <size_t k, class Data>
	void read(ifstream& in, multi_array<Data,k>& m) {
		range_traverser<k> t(reinterpret_cast<const array<int,k>&>(*m.shape()));
		do {
			in.read(reinterpret_cast<char *>(&m(*t)),sizeof(Data));
		} while(++t);
	}

	//! Writes a multi array to a stream
	template <size_t k, class Data>
	void write(ofstream& out, multi_array<Data,k>& m) {
		range_traverser<k> t(reinterpret_cast<const array<int,k>&>(*m.shape()));
		do {
			out.write(reinterpret_cast<char *>(&m(*t)),sizeof(Data));
		} while(++t);
	}

	//! Loads a multi array from a file
	template <class Data, size_t k>
	multi_array<Data,k> * load(char * filename) {
		ifstream in(filename,ios::binary|ios::in);
		array<int,k> a;
		read(in,a);
		multi_array<Data,k> * mp=new multi_array<Data,k>(a);
		read(in,*mp);
		return mp;
	}

	//! Saves a multi array to a file
	template <class Data, size_t k>
	void save(char * filename, multi_array<Data,k>& m) {
		ofstream out(filename,ios::binary|ios::out);
		array<int,k> a(reinterpret_cast<const array<int,k>&>(*m.shape()));
		write(out,a);
		write(out,m);
	}

	/*! @} */

}

#endif //VGTL_ARRAY_IO_HPP
