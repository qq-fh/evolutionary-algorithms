#ifndef RANDOM_H
#define RANDOM_H



#define randomInt(min, max)  ( min+rand()%(max-min) )   //随机整数
#define conv(x) ( (int)(1/x) )
#define randomDouble(min, max, precision) ( min+ (double)(rand()%conv(precision))/conv(precision)*(max-min)) //随机浮点数




#endif
