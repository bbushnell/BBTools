package bin;

import parse.LineParser1;

/**
 * Parses and manages Genome Taxonomy Database (GTDB) taxonomic lineage information.
 * Extracts taxonomic hierarchy from tab-delimited input with semicolon-separated classification strings.
 * Supports retrieval of complete taxonomic paths at different hierarchical levels.
 * @author Brian Bushnell
 */
public class GTDBLine {

	/**
	 * Constructs a GTDBLine from a tab-delimited byte array.
	 * Parses the line using tab and semicolon delimiters to extract taxonomic information.
	 * @param line Tab-delimited byte array containing name and taxonomic classification
	 */
	public GTDBLine(byte[] line) {
		//TODO: Possible bug [bin/GTDBLine#002] - dead ctor (no callers); passes an unset lpsemi -> the 2-arg ctor's lpsemi.parseString(1) would NPE/AIOOBE if ever called. Latent LOW.
		this(new LineParser1('\t').set(line), new LineParser1(';'));
	}
	
	/**
	 * Constructs a GTDBLine using pre-configured line parsers.
	 * Extracts name from tab parser and taxonomic hierarchy from semicolon parser.
	 * Parses GTDB format classification strings with prefixes (d__, p__, c__, o__, f__, g__, s__).
	 *
	 * @param lptab Tab-delimited line parser for extracting the sequence name
	 * @param lpsemi Semicolon-delimited parser for extracting taxonomic classification
	 */
	public GTDBLine(LineParser1 lptab, LineParser1 lpsemi) {
		//TODO: Possible bug [bin/GTDBLine#001] - javadoc says "pre-configured" parsers, but sole caller GradeBins.loadGTDBFile:1263 never set()s lpsemi -> parseString(1) NPEs. Also classification should be the FULL lineage = lptab.parseString(1) (tab field 1, fed to new Lineage()), not semi-term 1. MEDIUM, gated on gtdb= reachability; escalated. See report.
		name=lptab.parseString(0);//tab field 0 = user_genome name; lptab IS set by caller -> OK
		classification=lpsemi.parseString(1);
		if(classification==null) {return;}//dead guard: parseString throws on a null line, never returns null
		
		//d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Pacearchaeales;f__GW2011-AR1;g__MWBV01;s__
		//Dispatch on first char (d/p/c/o/f/g/s - all distinct GTDB prefixes); asserts one term per rank; any other first char silently skipped. charAt(0) safe iff every term is non-empty (GTDB X__ prefix always present; a trailing ';' empty term would crash, but that is format-impossible).
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
	
	/**
	 * Returns the complete taxonomic path up to the specified level.
	 * Level 0 returns domain only, higher levels return progressively longer paths.
	 *
	 * @param level Taxonomic level (0=domain, 1=phylum, 2=class, 3=order, 4=family, 5=genus, 6=species)
	 * @return Complete taxonomic path from domain to the specified level
	 * @throws RuntimeException if level is outside valid range 0-6
	 */
	//level 0 = raw domain term; 1-6 = cumulative path via phylum()..species(); out-of-range throws (crash-loud, correct). Requires domain..species populated -> only true if the term loop ran (lpsemi was set; see #001).
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

	//obvious: phylum()..species() each prepend the previous to build a cumulative domain;...;rank path. A null rank (term loop never ran) yields a literal "null;" segment rather than a crash - another #001 symptom.
	/** Returns the taxonomic path from domain to phylum level */
	public String phylum() {return domain+";"+phylum;}
	/** Returns the taxonomic path from domain to class level */
	public String classname() {return phylum()+";"+classname;}
	/** Returns the taxonomic path from domain to order level */
	public String order() {return classname()+";"+order;}
	/** Returns the taxonomic path from domain to family level */
	public String family() {return order()+";"+family;}
	/** Returns the taxonomic path from domain to genus level */
	public String genus() {return family()+";"+genus;}
	/** Returns the complete taxonomic path from domain to species level */
	public String species() {return genus()+";"+species;}
	
	/** Name or identifier of the genomic sequence */
	String name;
	
	/** Domain taxonomic level (prefix d__) */
	String domain;
	/** Phylum taxonomic level (prefix p__) */
	String phylum;
	/** Class taxonomic level (prefix c__) */
	String classname;
	/** Order taxonomic level (prefix o__) */
	String order;
	/** Family taxonomic level (prefix f__) */
	String family;
	/** Genus taxonomic level (prefix g__) */
	String genus;
	/** Species taxonomic level (prefix s__) */
	String species;
	/** Complete semicolon-delimited taxonomic classification string */
	String classification;
	
}
