#ifndef _CYCLEQUEUE_H
#define _CYCLEQUEUE_H
#include <iostream>

template <class T>
class cycleQueue
{
    private: 
        int m_size;
        int m_front;
        int m_rear;
        T*  m_data;
    public:
        cycleQueue(int size)
            :m_size(size),
            m_front(0),
            m_rear(0)
        {   
            m_data = new T[size];
        }   

        ~cycleQueue()
        {   
            delete [] m_data;
        }   

        bool isEmpty()
        {   
            return m_front == m_rear;
        }   

        bool isFull() 
        {   
            return m_front == (m_rear + 1)%m_size;
        }   

        void push(T ele)throw(std::bad_exception)
        {
            if(isFull()) {
                throw std::bad_exception();
            }
            m_data[m_rear] = ele;
            m_rear = (m_rear + 1)%m_size;
        }

        T pop() throw(std::bad_exception)
        {
            if(isEmpty())
            {
                throw std::bad_exception();
            }
            T tmp = m_data[m_front];
            m_front = (m_front + 1)%m_size;
            return tmp;
        }
};

/*
int main()
{
    cycleQueue<int> q(5);
    q.push(1);
    q.push(2);
    q.push(3);
    q.push(4);
    for (int i = 0; i < 4 ; i++){
        cout << q.pop() << endl;
		cout<<"-----i: "<<i <<endl;
		}
    q.push(5);
    q.push(5);
    q.push(5);
	q.push(77);
    cout << q.pop() << endl;
    cout << q.pop() << endl;
    cout << q.pop() << endl;
    cout << q.pop() << endl;
  	cout<< "----"<<endl;
  return 0;
}*/


#endif
