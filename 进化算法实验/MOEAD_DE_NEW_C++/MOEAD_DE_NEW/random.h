#ifndef RANDOM_H
#define RANDOM_H



#define randomInt(min, max)  ( min+rand()%(max-min) )   //�������
#define conv(x) ( (int)(1/x) )
#define randomDouble(min, max, precision) ( min+ (double)(rand()%conv(precision))/conv(precision)*(max-min)) //���������




#endif
