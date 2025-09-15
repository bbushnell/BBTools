package json;

import shared.Tools;

public class JsonLiteral {
	
	public JsonLiteral(String s_){
		s=s_;
	}
	
	public JsonLiteral(double value, int decimals){
		s=Tools.format("%."+decimals+"f", value);
	}
	
	@Override
	public String toString(){return s;}
	
	private final String s;
	
}
