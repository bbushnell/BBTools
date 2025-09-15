package fun;

public class ParseDouble {
    // Lookup table for decimal multipliers
    private static final double[] DECIMAL_INV_MULT = new double[20];
    static {
        for(int i=0; i<DECIMAL_INV_MULT.length; i++) {
            DECIMAL_INV_MULT[i] = Math.pow(0.1, i);
        }
    }

    // For special value checking
    private static final boolean[] NUMERIC_MAP = new boolean[128];
    static {
        for(int i='0'; i<='9'; i++) {
            NUMERIC_MAP[i] = true;
        }
        NUMERIC_MAP['-'] = true;
        NUMERIC_MAP['.'] = true;
    }

    public static void main(String[] args) {
        if(args.length<2) {
            System.out.println("Usage: DoubleParseTest <number> <iterations>");
            return;
        }
        
        String valueStr = args[0];
        int iterations = Integer.parseInt(args[1]);
        
        System.out.println("Testing with value: " + valueStr);
        System.out.println("Iterations: " + iterations);
        
        // Warmup for fair JIT optimization
        System.out.println("Warming up...");
        for(int i=0; i<5; i++) {
            testCustomParser(valueStr, iterations/10);
            testJavaParser(valueStr, iterations/10);
        }
        
        // Test custom parser
        long startTime = System.nanoTime();
        double customSum = testCustomParser(valueStr, iterations);
        long customTime = System.nanoTime() - startTime;
        
        // Test custom parser 2
        startTime = System.nanoTime();
        double customSum2 = testCustomParser2(valueStr, iterations);
        long customTime2 = System.nanoTime() - startTime;
        
        // Test Java parser
        startTime = System.nanoTime();
        double javaSum = testJavaParser(valueStr, iterations);
        long javaTime = System.nanoTime() - startTime;
        
        // Print results
        System.out.println("\nResults:");
        System.out.println("Custom parser time:  \t" + (customTime/1_000_000.0) + " ms");
        System.out.println("Custom parser time2: \t" + (customTime2/1_000_000.0) + " ms");
        System.out.println("Java parser time:    \t" + (javaTime/1_000_000.0) + " ms");
        System.out.println("Speed ratio:         \t" + (double)javaTime/customTime + "x");
        
        System.out.println("\nValidation sums:");
        System.out.println("Custom parser sum: " + customSum);
        System.out.println("Java parser sum: " + javaSum);
        System.out.println("Difference: " + Math.abs(customSum - javaSum));
    }
    
    private static double testCustomParser(String valueStr, int iterations) {
        byte[] bytes = valueStr.getBytes();
        double sum = 0;
        
        for(int i=0; i<iterations; i++) {
            sum += parseDouble(bytes, 0, bytes.length);
        }
        
        return sum;
    }
    
    private static double testCustomParser2(String valueStr, int iterations) {
        byte[] bytes = valueStr.getBytes();
        double sum = 0;
        
        for(int i=0; i<iterations; i++) {
            sum += parseDouble2(bytes, 0, bytes.length);
        }
        
        return sum;
    }
    
    private static double testJavaParser(String valueStr, int iterations) {
        double sum = 0;
        
        for(int i=0; i<iterations; i++) {
            sum += Double.parseDouble(valueStr);
        }
        
        return sum;
    }
    
    // Custom parseDouble implementation
    public static double parseDouble(final byte[] array, final int a0, final int b) {
        // Don't check for special values
        
        int a=a0;
        long upper=0;
        final byte z='0';
        long mult=1;
        if(array[a]=='-') {mult=-1; a++;}
        
        // Parse integer part
        for(; a<b; a++) {
            final byte c=array[a];
            if(c=='.') {break;}
            if(c=='e' || c=='E') {
                // Found exponent - use Java parser
                return Double.parseDouble(new String(array, a0, b-a0));
            }
            final int x=(c-z);
            upper=(upper*10)+x;
        }
        
        // Parse decimal part
        long lower=0;
        int places=0;
        for(a++; a<b; a++) {
            final byte c=array[a];
            if(c=='e' || c=='E') {
                // Found exponent - use Java parser
                return Double.parseDouble(new String(array, a0, b-a0));
            }
            final int x=(c-z);
            lower=(lower*10)+x;
            places++;
        }
        
        double d=mult*(upper+lower*(places<DECIMAL_INV_MULT.length?
                DECIMAL_INV_MULT[places]:Math.pow(0.1, places)));
        return d;
    }
    
    // Custom parseDouble implementation
    public static double parseDouble2(final byte[] array, final int a0, final int b) {
        // Check for special values
        if(b-a0>1 && b-a0<5) {
            final byte x=array[a0];
            if(!NUMERIC_MAP[x&0x7F]) {
                if(x=='N') {return Double.NaN;}
                if(x=='I') {return Double.POSITIVE_INFINITY;}
            }
            if(x=='-' && array[a0+1]=='I') {
                return Double.NEGATIVE_INFINITY;
            }
        }
        
        int a=a0;
        long upper=0;
        final byte z='0';
        long mult=1;
        if(array[a]=='-') {mult=-1; a++;}
        
        // Parse integer part
        for(; a<b; a++) {
            final byte c=array[a];
            if(c=='.') {break;}
            if(c=='e' || c=='E') {
                // Found exponent - use Java parser
                return Double.parseDouble(new String(array, a0, b-a0));
            }
            final int x=(c-z);
            upper=(upper*10)+x;
        }
        
        // Parse decimal part
        long lower=0;
        int places=0;
        for(a++; a<b; a++) {
            final byte c=array[a];
            if(c=='e' || c=='E') {
                // Found exponent - use Java parser
                return Double.parseDouble(new String(array, a0, b-a0));
            }
            final int x=(c-z);
            lower=(lower*10)+x;
            places++;
        }
        
        double d=mult*(upper+lower*(places<DECIMAL_INV_MULT.length?
                DECIMAL_INV_MULT[places]:Math.pow(0.1, places)));
        return d;
    }
}