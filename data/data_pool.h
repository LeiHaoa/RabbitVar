#ifndef DATA_POOL_H
#define DATA_POOL_H

#include <vector>
#include "../Variation.h"
class dataPool{
public:
	vector<Variation*> _data;
	int _pointer;
	int _size;

	dataPool(int init_size){
		cout << "init a variation with size: " << init_size << endl;
		//_data.resize(init_size);
		_data.reserve(init_size);
		_size = init_size;
		for(int i = 0; i < init_size; i++){
			_data.push_back(new Variation());
		}
	}
	void reset(){
		_pointer = 0;
	}
	Variation* get_variation(){
		if(_pointer == (_size - 1)){
			//_data.push_back(Variation());
			_data.reserve(_size + _size/2);
			for(int i = _size; i < _data.capacity(); i++){
				_data.push_back(new Variation());
			}
			_size = _size + _size/2;
			cerr << "resize data buf to: " << _size << "capacity: " << _data.capacity() << endl;
		}
		//cout << "get a variation: " << _pointer << " data: " << _data[_pointer].varsCountOnReverse << endl;
		Variation* var = _data[_pointer];
		//if(_pointer > 88835){
		//	//var->incDir(true) ;
		//	//cout << _pointer <<  " - after assigment:  " << var->varsCountOnReverse << endl;
		//}
		_pointer++;
		return var;
	}
};

#endif
