package bin;

import shared.LineParser1;

public class GTDBLine {

	public GTDBLine(byte[] line) {
		this(new LineParser1('\t').set(line), new LineParser1(';'));
	}
	
	public GTDBLine(LineParser1 lptab, LineParser1 lpsemi) {
		name=lptab.parseString(0);
		classification=lpsemi.parseString(1);
		if(classification==null) {return;}
		
		//d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__GW2011-AR1;g__MWBV01;s__
		for(int i=0; i<lpsemi.terms(); i++) {
			String term=lpsemi.parseString(i);
			char c=term.charAt(0);
			if(c=='d') {
				assert(domain==null) : domain+" -> "+term;
				domain=term;
			}else if(c=='p') {
				assert(phylum==null) : phylum+" -> "+term;
				phylum=term;
			}else if(c=='c') {
				assert(classname==null) : classname+" -> "+term;
				classname=term;
			}else if(c=='o') {
				assert(order==null) : order+" -> "+term;
				order=term;
			}else if(c=='f') {
				assert(family==null) : family+" -> "+term;
				family=term;
			}else if(c=='g') {
				assert(genus==null) : genus+" -> "+term;
				genus=term;
			}else if(c=='s') {
				assert(species==null) : species+" -> "+term;
				species=term;
			}
		}
	}
	
	public String getLevel(int level) {
		if(level==0) {return domain;}
		else if(level==1) {return phylum();}
		else if(level==2) {return classname();}
		else if(level==3) {return order();}
		else if(level==4) {return family();}
		else if(level==5) {return genus();}
		else if(level==6) {return species();}
		throw new RuntimeException("Bad level "+level);
	}

	public String phylum() {return domain+";"+phylum;}
	public String classname() {return phylum()+";"+classname;}
	public String order() {return classname()+";"+order;}
	public String family() {return order()+";"+family;}
	public String genus() {return family()+";"+genus;}
	public String species() {return genus()+";"+species;}
	
	String name;
	
	String domain;
	String phylum;
	String classname;
	String order;
	String family;
	String genus;
	String species;
	String classification;
	
}
