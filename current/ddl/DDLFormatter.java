package ddl;

import structures.ByteBuilder;

/**
 * Instance-based output formatter for DDL comparison results.
 * Controls column selection, output format, and result rendering.
 * All mutable state is per-instance — safe for multi-use, server mode,
 * and concurrent callers with different settings.
 *
 * Owns its own {@link #parse} method for CLI argument delegation:
 * the caller's parse loop calls {@code if(formatter.parse(arg, a, b)){}} first,
 * falling through to its own handling only if the formatter doesn't claim the arg.
 *
 * @author Noire, Brian Bushnell
 * @date May 17, 2026
 */
public class DDLFormatter implements Cloneable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DDLFormatter(){}

	@Override
	public DDLFormatter clone(){
		try{
			return (DDLFormatter)super.clone();
		}catch(CloneNotSupportedException e){
			throw new RuntimeException(e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses a single CLI argument.
	 * @param arg Original argument string
	 * @param a Lowercase parameter name (leading hyphens stripped)
	 * @param b Value string (may be null for boolean flags)
	 * @return true if this argument was recognized and handled */
	public boolean parse(String arg, String a, String b){
		if(a.equals("format")){
			if(b==null || b.equalsIgnoreCase("tab") || b.equalsIgnoreCase("tabular")){
				format=FORMAT_TAB;
			}else if(b.equalsIgnoreCase("json")){
				format=FORMAT_JSON;
			}else if(b.equalsIgnoreCase("oneline") || b.equalsIgnoreCase("3col")
					|| b.equalsIgnoreCase("query_ref_ani")){
				format=FORMAT_QUERY_REF_ANI;
			}else{
				format=Integer.parseInt(b);
			}
		}
		else if(a.equals("json")){
			if(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
				format=FORMAT_JSON;
			}
		}
		else if(a.equals("printani") || a.equals("ani")){printANI=parseBool(b);}
		else if(a.equals("printwkid") || a.equals("wkid")){printWKID=parseBool(b);}
		else if(a.equals("printcompleteness") || a.equals("completeness")){printCompleteness=parseBool(b);}
		else if(a.equals("printcontainment") || a.equals("containment")){printContainment=parseBool(b);}
		else if(a.equals("printmatches") || a.equals("matches")){printMatches=parseBool(b);}
		else if(a.equals("printbases") || a.equals("bases")){printBases=parseBool(b);}
		else if(a.equals("printtaxid") || a.equals("taxid") || a.equals("tid")){printTaxID=parseBool(b);}
		else if(a.equals("printname") || a.equals("name")){printName=parseBool(b);}
		else if(a.equals("printcardinality") || a.equals("cardinality") || a.equals("card")){printCardinality=parseBool(b);}
		else if(a.equals("printqueryname") || a.equals("queryname") || a.equals("qname")){printQueryName=parseBool(b);}
		else if(a.equals("printssu") || a.equals("ssuani")){printSSU=parseBool(b);}
		else if(a.equals("printall")){
			boolean x=parseBool(b);
			printANI=x; printWKID=x; printCompleteness=x; printContainment=x;
			printMatches=x; printBases=x; printTaxID=x; printName=x;
			printCardinality=x; printQueryName=x; printSSU=x;
		}
		else if(a.equals("noheader")){noHeader=parseBool(b);}
		else{return false;}
		return true;
	}

	private static boolean parseBool(String b){
		return b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Format Dispatch        ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats a header line into the ByteBuilder.
	 * Call once before formatting results (unless noHeader is set). */
	public void header(ByteBuilder bb){
		if(noHeader){return;}
		switch(format){
			case FORMAT_TAB: headerTab(bb); break;
			case FORMAT_JSON: break;
			case FORMAT_QUERY_REF_ANI: headerQueryRefAni(bb); break;
			default: headerTab(bb); break;
		}
	}

	/** Formats a single comparison result into the ByteBuilder.
	 * @param c The comparison to format
	 * @param bb Destination ByteBuilder */
	public void format(DDLComparison c, ByteBuilder bb){
		switch(format){
			case FORMAT_TAB: formatTab(c, bb); break;
			case FORMAT_JSON: formatJson(c, bb); break;
			case FORMAT_QUERY_REF_ANI: formatQueryRefAni(c, bb); break;
			default: formatTab(c, bb); break;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Tabular Format         ----------------*/
	/*--------------------------------------------------------------*/

	private void headerTab(ByteBuilder bb){
		boolean first=true;
		if(printANI){if(!first){bb.tab();} bb.append("ANI"); first=false;}
		if(printWKID){if(!first){bb.tab();} bb.append("WKID"); first=false;}
		if(printCompleteness){if(!first){bb.tab();} bb.append("Complt"); first=false;}
		if(printContainment){if(!first){bb.tab();} bb.append("Cont"); first=false;}
		if(printMatches){if(!first){bb.tab();} bb.append("Matches"); first=false;}
		if(printBases){if(!first){bb.tab();} bb.append("Bases"); first=false;}
		if(printCardinality){if(!first){bb.tab();} bb.append("Card"); first=false;}
		if(printSSU){if(!first){bb.tab();} bb.append("SSU"); first=false;}
		if(printTaxID){if(!first){bb.tab();} bb.append("TID"); first=false;}
		if(printQueryName){if(!first){bb.tab();} bb.append("Query"); first=false;}
		if(printName){if(!first){bb.tab();} bb.append("Name"); first=false;}
		bb.nl();
	}

	private void formatTab(DDLComparison c, ByteBuilder bb){
		boolean first=true;
		if(printANI){if(!first){bb.tab();} appendFloat(bb, c.ani, 4); first=false;}
		if(printWKID){if(!first){bb.tab();} appendMetric(bb, c.wkid, 4); first=false;}
		if(printCompleteness){if(!first){bb.tab();} appendMetric(bb, c.completenessAB, 4); first=false;}
		if(printContainment){if(!first){bb.tab();} appendMetric(bb, c.containmentAB, 4); first=false;}
		if(printMatches){if(!first){bb.tab();} bb.append(c.equal); first=false;}
		if(printBases){if(!first){bb.tab();} bb.append(c.refRecord!=null ? c.refRecord.bases : 0); first=false;}
		if(printCardinality){if(!first){bb.tab();} bb.append(c.refRecord!=null ? c.refRecord.cardinality : -1); first=false;}
		if(printSSU){if(!first){bb.tab();} appendMetric(bb, c.ssuIdentity, 4); first=false;}
		if(printTaxID){if(!first){bb.tab();} bb.append(c.refRecord!=null ? c.refRecord.taxID : -1); first=false;}
		if(printQueryName){if(!first){bb.tab();} bb.append(qname(c)); first=false;}
		if(printName){if(!first){bb.tab();} bb.append(rname(c)); first=false;}
		bb.nl();
	}

	/*--------------------------------------------------------------*/
	/*----------------        JSON Format           ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts a JSON array. Call before formatting individual results. */
	public void jsonStart(ByteBuilder bb){
		bb.append('[').nl();
		jsonFirst=true;
	}

	/** Ends a JSON array. Call after all results are formatted. */
	public void jsonEnd(ByteBuilder bb){
		bb.nl().append(']').nl();
	}

	private void formatJson(DDLComparison c, ByteBuilder bb){
		if(!jsonFirst){bb.append(',').nl();}
		jsonFirst=false;
		bb.append("  {");
		boolean f=true;
		if(printANI){if(!f){bb.append(',');} bb.append("\"ANI\":"); appendJsonFloat(bb, c.ani, 6); f=false;}
		if(printWKID && c.wkid>=0){if(!f){bb.append(',');} bb.append("\"WKID\":"); appendJsonFloat(bb, c.wkid, 6); f=false;}
		if(printCompleteness && c.completenessAB>=0){if(!f){bb.append(',');} bb.append("\"Completeness\":"); appendJsonFloat(bb, c.completenessAB, 6); f=false;}
		if(printContainment && c.containmentAB>=0){if(!f){bb.append(',');} bb.append("\"Containment\":"); appendJsonFloat(bb, c.containmentAB, 6); f=false;}
		if(printMatches){if(!f){bb.append(',');} bb.append("\"Matches\":").append(c.equal); f=false;}
		if(printBases && c.refRecord!=null){if(!f){bb.append(',');} bb.append("\"Bases\":").append(c.refRecord.bases); f=false;}
		if(printCardinality && c.refRecord!=null){if(!f){bb.append(',');} bb.append("\"Cardinality\":").append(c.refRecord.cardinality); f=false;}
		if(printSSU && c.ssuIdentity>=0){if(!f){bb.append(',');} bb.append("\"SSU\":"); appendJsonFloat(bb, c.ssuIdentity, 4); f=false;}
		if(printTaxID && c.refRecord!=null){if(!f){bb.append(',');} bb.append("\"TaxID\":").append(c.refRecord.taxID); f=false;}
		String qn=qname(c); if(printQueryName && !qn.equals("-")){if(!f){bb.append(',');} bb.append("\"Query\":\""); escapeJson(bb, qn); bb.append('"'); f=false;}
		String rn=rname(c); if(printName && !rn.equals("-")){if(!f){bb.append(',');} bb.append("\"Name\":\""); escapeJson(bb, rn); bb.append('"'); f=false;}
		bb.append('}');
	}

	private static void escapeJson(ByteBuilder bb, String s){
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='"'){bb.append('\\').append('"');}
			else if(c=='\\'){bb.append('\\').append('\\');}
			else if(c=='\n'){bb.append('\\').append('n');}
			else if(c=='\t'){bb.append('\\').append('t');}
			else{bb.append(c);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Query-Ref-ANI Format      ----------------*/
	/*--------------------------------------------------------------*/

	private void headerQueryRefAni(ByteBuilder bb){
		bb.append("Query\tRef\tANI");
		if(printWKID){bb.append("\tWKID");}
		if(printCompleteness){bb.append("\tComplt");}
		if(printMatches){bb.append("\tMatches");}
		if(printSSU){bb.append("\tSSU");}
		if(printTaxID){bb.append("\tTID");}
		bb.nl();
	}

	private void formatQueryRefAni(DDLComparison c, ByteBuilder bb){
		bb.append(qname(c));
		bb.tab();
		bb.append(rname(c));
		bb.tab();
		appendFloat(bb, c.ani, 4);
		if(printWKID){bb.tab(); appendMetric(bb, c.wkid, 4);}
		if(printCompleteness){bb.tab(); appendMetric(bb, c.completenessAB, 4);}
		if(printMatches){bb.tab(); bb.append(c.equal);}
		if(printSSU){bb.tab(); appendMetric(bb, c.ssuIdentity, 4);}
		if(printTaxID){bb.tab(); bb.append(c.refRecord!=null ? c.refRecord.taxID : -1);}
		bb.nl();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	/** Query name from record, or "-" if unavailable. */
	private static String qname(DDLComparison c){
		return c.queryRecord!=null && c.queryRecord.name!=null ? c.queryRecord.name : "-";
	}

	/** Ref name from record, or "-" if unavailable. */
	private static String rname(DDLComparison c){
		return c.refRecord!=null && c.refRecord.name!=null ? c.refRecord.name : "-";
	}

	/** Appends a float with fixed decimal places. */
	private static void appendFloat(ByteBuilder bb, float v, int decimals){
		bb.appendSlow(v, decimals);
	}

	/** Appends a metric that may be unavailable (-1). Prints "-" if unavailable. */
	private static void appendMetric(ByteBuilder bb, float v, int decimals){
		if(v<0){bb.append('-');}
		else{bb.appendSlow(v, decimals);}
	}

	/** Appends a float for JSON (no padding). */
	private static void appendJsonFloat(ByteBuilder bb, float v, int decimals){
		bb.appendSlow(v, decimals);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/* Format selection */
	public int format=FORMAT_TAB;

	/* Column visibility — all instance fields */
	public boolean printANI=true;
	public boolean printWKID=true;
	public boolean printCompleteness=true;
	public boolean printContainment=false;
	public boolean printMatches=true;
	public boolean printBases=true;
	public boolean printTaxID=true;
	public boolean printName=true;
	public boolean printQueryName=false;
	public boolean printCardinality=false;
	public boolean printSSU=false;

	/* Header control */
	public boolean noHeader=false;

	/* JSON state */
	private boolean jsonFirst=true;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int FORMAT_TAB=0;
	public static final int FORMAT_JSON=1;
	public static final int FORMAT_QUERY_REF_ANI=2;
}
