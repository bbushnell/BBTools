package cardinality;
public class MergeTest {
    public static void main(String[] args) {
        // ErtlULL
        ErtlULL ea=new ErtlULL(2048, 31, -1, 0);
        ErtlULL eb=new ErtlULL(2048, 31, -1, 0);
        ErtlULL ec=new ErtlULL(2048, 31, -1, 0);
        for(long i=0; i<100000; i++){
            if(i%2==0){ea.add(i);}else{eb.add(i);}
            ec.add(i);
        }
        ea.add(eb);
        System.out.println("ErtlULL: merged="+ea.cardinality()+" expected="+ec.cardinality());

        // UDLL6
        UltraDynamicLogLog6 ua=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6 ub=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6 uc=new UltraDynamicLogLog6(2048, 31, -1, 0);
        for(long i=0; i<100000; i++){
            if(i%2==0){ua.add(i);}else{ub.add(i);}
            uc.add(i);
        }
        ua.add(ub);
        System.out.println("UDLL6:   merged="+ua.cardinality()+" expected="+uc.cardinality());

        // Higher cardinality
        ErtlULL ha=new ErtlULL(2048, 31, -1, 0);
        ErtlULL hb=new ErtlULL(2048, 31, -1, 0);
        ErtlULL hc=new ErtlULL(2048, 31, -1, 0);
        for(long i=0; i<10000000; i++){
            if(i%2==0){ha.add(i);}else{hb.add(i);}
            hc.add(i);
        }
        ha.add(hb);
        System.out.println("ErtlULL 10M: merged="+ha.cardinality()+" expected="+hc.cardinality());

        UltraDynamicLogLog6 da=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6 db=new UltraDynamicLogLog6(2048, 31, -1, 0);
        UltraDynamicLogLog6 dc=new UltraDynamicLogLog6(2048, 31, -1, 0);
        for(long i=0; i<10000000; i++){
            if(i%2==0){da.add(i);}else{db.add(i);}
            dc.add(i);
        }
        da.add(db);
        System.out.println("UDLL6 10M:   merged="+da.cardinality()+" expected="+dc.cardinality());
    }
}
