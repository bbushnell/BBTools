package test.parse;

import parse.Parse;

/** Timing fixture for byte/String ASCII-insensitive comparison candidates. */
public final class ParseByteStringCaseTimingTest{

	/** Runs one implementation and corpus shape with no allocation in the timed loop.
	 * @param args implementation=production|fixed|guarded|table|shorttable|hybrid,
	 * shape=exact|mixed|early|late, length=32|256, plus optional rounds/warmup */
	public static void main(String[] args){
		String implementation="guarded", shape="mixed";
		int length=32;
		long rounds=20000, warmup=2000;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Bad argument "+arg);}
			final String key=split[0], value=split[1];
			if(key.equals("implementation")){implementation=value;}
			else if(key.equals("shape")){shape=value;}
			else if(key.equals("length")){length=Integer.parseInt(value);}
			else if(key.equals("rounds")){rounds=Long.parseLong(value);}
			else if(key.equals("warmup")){warmup=Long.parseLong(value);}
			else{throw new IllegalArgumentException("Unknown argument "+arg);}
		}
		require(length==32 || length==256, "length must be 32 or 256");
		require(rounds>=1 && warmup>=0, "round and warmup counts must be nonnegative");
		if(implementation.equals("production") && !shape.equals("exact")){
			throw new IllegalArgumentException(
					"production timing is only a fair baseline for exact-case inputs");
		}
		final Corpus corpus=new Corpus(shape, length, 1024);
		run(implementation, corpus, warmup);
		final long start=System.nanoTime();
		final long sum=run(implementation, corpus, rounds);
		final long elapsed=System.nanoTime()-start;
		final long expected=corpus.expectedChecksum*rounds;
		require(sum==expected, "checksum "+sum+" != "+expected);
		System.out.println("BYTE_STRING_CASE_TIMING_RESULT implementation="+implementation+
				" shape="+shape+" length="+length+" rounds="+rounds+" warmup="+warmup+
				" pairs="+(rounds*corpus.size())+" nanos="+elapsed+" checksum="+sum);
	}

	/** Dispatches once to a loop dedicated to the selected implementation. */
	private static long run(String implementation, Corpus corpus, long rounds){
		if(implementation.equals("production") || implementation.equals("fixed")){
			return loopProduction(corpus, rounds);
		}
		if(implementation.equals("guarded")){return loopGuarded(corpus, rounds);}
		if(implementation.equals("table")){return loopTable(corpus, rounds);}
		if(implementation.equals("shorttable")){return loopShortTable(corpus, rounds);}
		if(implementation.equals("hybrid")){return loopHybrid(corpus, rounds);}
		throw new IllegalArgumentException("Unknown implementation "+implementation);
	}

	/** Current production byte/String comparison. */
	private static long loopProduction(Corpus corpus, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(int i=0; i<corpus.size(); i++){
				sum+=Parse.equalsIgnoreCase(corpus.bytes[i], corpus.strings[i]) ? i+1 : -i-1;
			}
		}
		return sum;
	}

	/** Exact-match short-circuit plus guarded OR32 letter folding. */
	private static long loopGuarded(Corpus corpus, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(int i=0; i<corpus.size(); i++){
				sum+=equalsGuarded(corpus.bytes[i], corpus.strings[i]) ? i+1 : -i-1;
			}
		}
		return sum;
	}

	/** Unconditional ASCII fold-table lookup after ASCII bounds checks. */
	private static long loopTable(Corpus corpus, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(int i=0; i<corpus.size(); i++){
				sum+=equalsTable(corpus.bytes[i], corpus.strings[i]) ? i+1 : -i-1;
			}
		}
		return sum;
	}

	/** Exact-match short-circuit plus ASCII fold-table lookup. */
	private static long loopShortTable(Corpus corpus, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(int i=0; i<corpus.size(); i++){
				sum+=equalsShortTable(corpus.bytes[i], corpus.strings[i]) ? i+1 : -i-1;
			}
		}
		return sum;
	}

	/** Exact-prefix scan followed by ASCII fold-table lookup. */
	private static long loopHybrid(Corpus corpus, long rounds){
		long sum=0;
		for(long r=0; r<rounds; r++){
			for(int i=0; i<corpus.size(); i++){
				sum+=equalsHybrid(corpus.bytes[i], corpus.strings[i]) ? i+1 : -i-1;
			}
		}
		return sum;
	}

	/** Candidate with exact-match short-circuit and guarded OR32 folding. */
	static boolean equalsGuarded(byte[] a, String b){
		if(a==null || b==null || a.length!=b.length()){return false;}
		for(int i=0; i<a.length; i++){
			final int x=a[i];
			final char c=b.charAt(i);
			if(x==c){continue;}
			final int folded=x|32;
			if(folded!=(c|32) || folded<'a' || folded>'z'){return false;}
		}
		return true;
	}

	/** Candidate with unconditional ASCII fold-table lookup. */
	static boolean equalsTable(byte[] a, String b){
		if(a==null || b==null || a.length!=b.length()){return false;}
		for(int i=0; i<a.length; i++){
			final int x=a[i];
			final char c=b.charAt(i);
			if(x<0 || c>127 || FOLD[x]!=FOLD[c]){return false;}
		}
		return true;
	}

	/** Candidate with exact-match short-circuit and ASCII fold-table lookup. */
	static boolean equalsShortTable(byte[] a, String b){
		if(a==null || b==null || a.length!=b.length()){return false;}
		for(int i=0; i<a.length; i++){
			final int x=a[i];
			final char c=b.charAt(i);
			if(x==c){continue;}
			if(x<0 || c>127 || FOLD[x]!=FOLD[c]){return false;}
		}
		return true;
	}

	/** Candidate with an exact-prefix scan and table lookup after the first mismatch. */
	static boolean equalsHybrid(byte[] a, String b){
		if(a==null || b==null || a.length!=b.length()){return false;}
		for(int i=0; i<a.length; i++){
			if(a[i]!=b.charAt(i)){return equalsFoldedFrom(a, b, i);}
		}
		return true;
	}

	/** Table-folds byte/String pairs from the first known mismatch. */
	private static boolean equalsFoldedFrom(byte[] a, String b, int from){
		for(int i=from; i<a.length; i++){
			final int x=a[i];
			final char c=b.charAt(i);
			if(x<0 || c>127 || FOLD[x]!=FOLD[c]){return false;}
		}
		return true;
	}

	/** Deterministic prebuilt byte/String inputs and checksum oracle. */
	private static final class Corpus{
		Corpus(String shape_, int length_, int size_){
			shape=shape_;
			bytes=new byte[size_][];
			strings=new String[size_];
			long checksum=0;
			for(int i=0; i<size_; i++){
				final byte[] array=new byte[length_];
				final char[] chars=new char[length_];
				for(int j=0; j<length_; j++){
					final int bits=mix(i*length_+j+0x9e3779b9);
					final byte b=base(bits);
					array[j]=b;
					chars[j]=paired(b, bits);
					if(shape.equals("exact")){chars[j]=(char)b;}
				}
				if(shape.equals("early")){chars[0]=(char)mismatch(array[0]);}
				else if(shape.equals("late")){chars[length_-1]=(char)mismatch(array[length_-1]);}
				else if(!shape.equals("exact") && !shape.equals("mixed")){
					throw new IllegalArgumentException("Unknown shape "+shape);
				}
				bytes[i]=array;
				strings[i]=new String(chars);
				checksum+=(shape.equals("exact") || shape.equals("mixed")) ? i+1 : -i-1;
			}
			expectedChecksum=checksum;
		}

		int size(){return bytes.length;}

		/** Shape used to build this corpus. */
		final String shape;
		/** Byte-side corpus. */
		final byte[][] bytes;
		/** String-side corpus. */
		final String[] strings;
		/** Expected checksum for one round. */
		final long expectedChecksum;
	}

	/** Varied deterministic ASCII letter/digit/punctuation value. */
	private static byte base(int bits){
		final int x=(bits>>>16)&63;
		if(x<26){return (byte)('A'+x);}
		if(x<52){return (byte)('a'+x-26);}
		if(x<62){return (byte)('0'+x-52);}
		return PUNCT[((bits>>>22)&0x7fffffff)%PUNCT.length];
	}

	/** Matching String-side char with deterministic letter case variation. */
	private static char paired(byte b, int bits){
		return isLetter(b) && (bits&1)==0 ? (char)(b^32) : (char)b;
	}

	/** ASCII mismatch that is not a case-folded match. */
	private static byte mismatch(byte b){
		return (b=='A' || b=='a') ? (byte)'B' : (byte)'A';
	}

	/** 32-bit deterministic integer mixer. */
	private static int mix(int x){
		x^=x>>>16;
		x*=0x7feb352d;
		x^=x>>>15;
		x*=0x846ca68b;
		x^=x>>>16;
		return x;
	}

	/** True for an ASCII letter. */
	private static boolean isLetter(byte b){
		return (b>='A' && b<='Z') || (b>='a' && b<='z');
	}

	/** Fails an expectation even when the JVM runs with -da. */
	private static void require(boolean condition, String reason){
		if(!condition){throw new AssertionError(reason);}
	}

	/** ASCII letter fold table. */
	private static final byte[] FOLD=makeFold();

	/** Punctuation characters available to the deterministic corpus. */
	private static final byte[] PUNCT={
		'!', '#', '$', '%', '&', '(', ')', '*', '+', ',', '-', '.', '/', ':', ';',
		'<', '=', '>', '?', '@', '[', ']', '^', '_', '`', '{', '|', '}', '~'
	};

	/** Creates the ASCII fold table used by the table candidates. */
	private static byte[] makeFold(){
		final byte[] table=new byte[128];
		for(int i=0; i<table.length; i++){table[i]=(byte)i;}
		for(int i='A'; i<='Z'; i++){table[i]=(byte)(i+32);}
		return table;
	}
}
