package fun;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.regex.Pattern;

import shared.Parse;

public class Foo2 {
	
	static final Pattern pipePattern=Pattern.compile("\\|");

	public static void main(String[] args) throws Exception{
		final boolean slow=args.length<2 ? true : "slow".equalsIgnoreCase(args[1]);
		final BufferedReader reader=new BufferedReader(new FileReader(args[0]));
		long sum=0, lines=0, chars=0;

//		ArrayList<Pair> list=new ArrayList<Pair>(10000000);
		for(String line=reader.readLine(); line!=null; line=reader.readLine()) {
			final long size=(slow ? processSlow(line) : processFast(line));
			if(size>=0) {
				sum+=size;
				lines++;
				chars+=line.length();
			}
		}

		reader.close();
		System.out.println("sum="+sum);
		System.out.println("lines="+lines);
		System.out.println("chars="+chars);
	}
	
	static long processSlow(String line) {
		String[] split=pipePattern.split(line);
		if(split[6].charAt(0)!='F') {return -1;}
		long size=Long.parseLong(split[3]);
		assert(size>=0) : line;
		return size;
	}

	static long processFast(String line) {
		final int delimiter='|';
		final int len=line.length();
		int a=0, b=0;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		//		long w=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		long size=Parse.parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term : '"+new String(line)+"'";
		if(line.charAt(a)!='F') {return -1;}
		b++;
		a=b;

		return size;
	}

}
