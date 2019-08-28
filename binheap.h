#ifndef BINHEAP_H
#define BINHEAP_H

/* Declarations */
template <class T>
class BinHeap;

/* Dependencies */
#include <iostream>
#include <math.h>

template <class T>
class BinHeap
{
public:
	/// Create a heap with a maximum number of elements.
	BinHeap(int max_size);

	/// Add a new element to the Heap and return the position.
	virtual int add(T x);

	/// Get the element with the minimum key.
	T get_min() const;

	/// Remove the element with the minimum key. Return this element.
	virtual T remove_min();

	/// Remove the element in position 'i'.
	virtual T remove(int i);

	/// Get the number of elements in the heap
	int size() const {return n;}

	/// Get the capacity of the heap;
	int max_size() const {return max_n;}

	/// Return true if the heap has no elements.
	bool is_empty() const {return n==0;}

	/// Return true if the number of elements is equal to the max size.
	bool is_full() const {return max_n == n;}

	/// Print the elements array.
	void print()const;

	/// Check if the heap is well formed.
	bool check() const;

	/// Check if the heap is well formed starting from node i.
	bool check(int i) const;

protected:
	T *elems;   /// array of heap elements
	int n;      /// number of elements
	int max_n;  /// maximum number of elements

	/// Move an element up in the heap until the heap shape property is fullfilled.
	int up_heap(int i);

	/// Move an element up down the heap until the heap shape property is fullfilled.
	int down_heap(int i);

	/// Computes the height of the Heap Tree.
	int height()const;

	/// Swap the elements in positions i and j.
	virtual void swap(int i, int j);

	/// This functions may be redefined in derived classes depending on T.
	virtual bool great_than(T elema, T elemb) const { return elema > elemb;}
        virtual bool less_than(T elema, T elemb) const { return elema < elemb;}
        virtual bool equal_to(T elema, T elemb) const { return elema == elemb;}
};

template <class T>
BinHeap<T>::BinHeap(int max_size)
{
	max_n = max_size;
	n = 0;
	if (max_size > 0) elems = new T[max_size];
	else elems = 0;
}

template <class T>
int BinHeap<T>::add(T x)
{
	if (is_full())
	{
		std::cout << "heap is full ... can not add a new element" << std::endl;
		return -1;
	}
	else
	{
                elems[n] = x; // insertar el elemento al final
		n++;
		return up_heap(n-1); // ponerlo en posicion
	}
}

template <class T>
T BinHeap<T>::remove_min()
{
	if (is_empty())
	{ 
		std::cout << "heap is empty ... can not remove min" << std::endl;
		throw 0;
	}
	else
	{
		swap(0, n-1);  // intercambiar el ultimo con el primero
		T x = elems[n-1]; // eliminar el elemento raiz que ahora esta en la ultima posicion
		n--;
                if (n > 0) down_heap(0); // poner la raiz en su lugar
		return x;
	}
}

template <class T>
T BinHeap<T>::remove(int i)
{
	if (i < 0 || i >= n)
	{
		std::cout << "trying to remove from a position out of range ..." << std::endl; 
		return 0;
	}

	T x;
	if (i == n-1) // estoy quitando el ultimo
	{
		x = elems[n-1];
		n--;
	}
	else
	{
        swap(i, n-1);
		x = elems[n-1]; // intercambiar la posicion 'i' con el ultimo
		n--;
		down_heap(i);
	}

	return x;
}

template <class T>
T BinHeap<T>::get_min() const
{
	if (is_empty())
	{
		std::cout << "heap is empty ... can not return min" << std::endl; 
		return 0;
	}
	return elems[0];
}


template <class T>
int BinHeap<T>::up_heap(int i)
{
	if (i < 0 || i >= n)
	{
		std::cout << "trying to make 'up_heap' from a position out of range ..." << std::endl; 
		throw 0;
		return -1;
	}
	
	while (i > 0 && less_than(elems[i], elems[(i-1)/2]))
	{
		swap((i-1)/2, i); // intercambiar hijo con padre
		i=(i-1)/2;        // moverme al padre
	}

	return i;
}

template <class T>
int BinHeap<T>::down_heap(int i)
{
	if (i < 0 || i >= n) 
	{
		std::cout << "trying to make 'down_heap' from a position out of range ..." << std::endl; 
		throw 0;
		return -1;
	}

	int h = height();
	int k = (int)pow((double)2, (int)h)-1;
	if (i >= k) return -1; // significa q es una hoja, nada que hacer
	
	while (i < k)
	{
		// si tiene dos hijos y si mi valor es mayor que los dos valores 
		// intercambiar con el menor de los dos
		if (2*i+2 < n && great_than(elems[i], elems[2*i+1]) && great_than(elems[i], elems[2*i+2]))
		{
			if (less_than(elems[2*i+1], elems[2*i+2])) // con el izquierdo
			{
                swap(2*i+1, i);
                i = 2*i+1;
			}
			else // con el derecho
			{
                swap(2*i+2, i);
                i = 2*i+2;
			}
		}
		else if (2*i+1 < n && great_than(elems[i], elems[2*i+1])) // si es mayor que el izquierdo
		{
			swap(2*i+1, i);
			i = 2*i+1;
		}
		else if (2*i+2 < n && great_than(elems[i], elems[2*i+2])) // si es mayor que el derecho
		{
			swap(2*i+2, i);
			i = 2*i+2;
		}
		else break;
	}
	return i;
}

template <class T>
int BinHeap<T>::height()const
{
	if (n == 0) return -1;
	
	int i = 0, h = 0;
	while (2*i+1 < n)
	{
		i = 2*i+1;
		h++;
	}

	return h;
}

template <class T>
void BinHeap<T>::swap(int i, int j)
{
	T temp = elems[i];
	elems[i] = elems[j];
	elems[j] = temp;
}

template <class T>
void BinHeap<T>::print()const
{
//	for (int i = 0; i < _size; i++)
//		std::cout << elems[i] << " ";
//	std::cout << std::endl;
}

template <class T>
bool BinHeap<T>::check() const
{
	return check(0);
}

template <class T>
bool BinHeap<T>::check(int i) const
{
	if (2*i+1 < n)
	{
		if (great_than(elems[i], elems[2*i+1])) return false;
		return check(2*i+1);
	}
	if (2*i+2 < n)
	{
		if (great_than(elems[i], elems[2*i+2])) return false;
		return check(2*i+2);
	}
	return true;
}

#endif
