package bin;

import json.JsonObject;
import shared.LineParserS2;
import structures.ByteBuilder;

public class SketchRecord extends BinObject {
	
	public SketchRecord(JsonObject hit) {
		setFrom(hit);
	}
	
	public void setFrom(JsonObject hit) {
		matches=hit.getLong("Matches").intValue();
		ani=hit.getDouble("ANI").floatValue();
		completeness=hit.getDouble("Complt").floatValue();
		contam=hit.getDouble("Contam").floatValue();
		taxid=hit.getLong("TaxID").intValue();
		taxName=hit.getString("taxName");
		genusTaxid=-1;//TODO, get from tree
	}
	
	public ByteBuilder toBytes() {
		ByteBuilder bb=new ByteBuilder();
		return appendTo(bb);
	}
	
	public ByteBuilder appendTo(ByteBuilder bb) {
		bb.append("ANI: ").append(ani, 2);
		bb.tab().append("Complt: ").append(completeness, 2);
		bb.tab().append("Contam: ").append(contam, 2);
		bb.tab().append("Matches: ").append(matches);
		bb.tab().append("TaxID: ").append(taxid);
		bb.tab().append("Name: ").append(shrink(taxName));
		return bb;
	}
	
	private static String shrink(String n) {
		if(n==null || n.length()<40) {return n;}
		ByteBuilder bb=new ByteBuilder();
		LineParserS2 lp=new LineParserS2(' ');
		lp.set(n);
		for(int i=0; i<4 && lp.hasMore(); i++) {
			if(bb.length>0) {bb.space();}
			bb.append(lp.parseString());
		}
		if(lp.hasMore()) {bb.append("...");}
		return (bb.length()<n.length() ? bb.toString() : n);
	}

	int matches=-1;
	float completeness=-1;
	float contam=-1;
	float ani=-1;
	int taxid=-1;
	int genusTaxid=-1;
	String taxName=null;
	
}
