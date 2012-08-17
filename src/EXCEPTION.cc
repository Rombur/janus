#include "EXCEPTION.hh"

void Check(bool condition,string message)
{
  if (condition==false)
    throw message;
}
