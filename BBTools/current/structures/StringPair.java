package structures;

public class StringPair {
	
	public StringPair(String a_, String b_){
		a=a_;
		b=b_;
	}
	
	@Override
	public String toString(){return "("+a+", "+b+")";}
	
	public String a;
	public String b;
	
}
