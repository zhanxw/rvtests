#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Argument.h"
#include "IO.h"
#include "tabix.h"
#include "PeopleSet.h"
#include "RangeList.h"
#include "Utils.h"
#include "VCFUtil.h"

#define COUNT 1000000

int main()
{
    const char a[] = "test.vcf.gz";
    VCFInputFile vin(a);

    for(int i=0; i<COUNT; i++)
    {
        vin.setRange("1", 0, 100);

        while(vin.readRecord())
        {
            VCFRecord &r = vin.getVCFRecord();
            VCFPeople &people = r.getPeople();
        }
    }
    return 0;
}
