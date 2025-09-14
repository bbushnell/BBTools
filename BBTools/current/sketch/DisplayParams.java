package sketch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map.Entry;

import json.JsonObject;
import shared.Colors;
import shared.Parse;
import shared.Tools;
import structures.ByteBuilder;
import tax.PrintTaxonomy;
import tax.TaxFilter;
import tax.TaxNode;
import tax.TaxTree;

public class DisplayParams implements Cloneable {
	
	@Override
	public DisplayParams clone(){
		try {
			DisplayParams copy=(DisplayParams)super.clone();
			if(taxFilterWhite!=null){
				copy.taxFilterWhite=taxFilterWhite.deepCopy();
			}
			if(taxFilterBlack!=null){
				copy.taxFilterBlack=taxFilterBlack.deepCopy();
			}
			copy.postParsed=false;
			return copy;
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	public DisplayParams parseDoubleHeader(String s){
		if(!s.startsWith("##")){return this;}
//		if(!s.startsWith("##")){return this.clone();}
		StringBuilder sb=new StringBuilder();
		for(int i=2; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='\n'){break;}
			sb.append(c);
		}
		return parseDoubleHeaderLine(sb.toString());
	}
	
	public DisplayParams parseDoubleHeaderLine(String line) {
		if(line.startsWith("##")){line=line.substring(2);}
		else{assert(!line.startsWith("#")) : line;}
		if(line.length()<1){return this;}
		
		DisplayParams params=this.clone();
		
		String[] args=line.split(" ");
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;} //Normally handled by PreParser, but not in this case.
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			boolean x=params.parse(arg, a, b);
//			assert(x) : "Unknown parameter "+arg+"\n"+line;
			if(!x){System.err.println("Warning: Unknown parameter "+arg);}
		}
		if(SketchObject.verbose2){System.err.println("Made it to post-parse.  taxFilterWhite="+params.taxFilterWhite);}
		params.postParse(true, true);
		if(SketchObject.verbose2){System.err.println("Passed post-parse.  taxFilterWhite="+params.taxFilterWhite);}
		
		return params;
	}
	
	public boolean parse(String arg, String a, String b){
	
		if(a.equals("chunk")){
			chunkNum=Integer.parseInt(b);
		}else if(a.equals("minhits")  || a.equals("hits")){
			minHits=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minwkid") || a.equalsIgnoreCase("wkid")){
			minWKID=Float.parseFloat(b);
			if(minWKID>1){minWKID/=100;}
			assert(minWKID<=1) : "minWKID should between 0 and 1";
		}else if(a.equalsIgnoreCase("minid") || a.equalsIgnoreCase("id") || a.equalsIgnoreCase("minani") || a.equalsIgnoreCase("ani")){
			minANI=Float.parseFloat(b);
			if(minANI>1){minANI/=100;}
			assert(minANI<=1) : "minANI should between 0 and 1";
			if(minANI>0){
				minWKID=(float)Tools.max(minWKID, Comparison.aniToWkid(minANI, 32));//Lowest possible minWKID for this ANI
			}
		}else if(a.equals("minbases")){
			minBases=Integer.parseInt(b);
		}else if(a.equals("minsizeratio")){
			minSizeRatio=Float.parseFloat(b);
//			assert(minSizeRatio>=0f && minSizeRatio<=1.0f) : "\nminSizeRatio must be between 0 and 1, inclusive.\n";
			if(minSizeRatio>1){minSizeRatio=1f/minSizeRatio;}
		}else if(a.equals("records") || a.equals("maxrecords") || a.equals("results")){
			maxRecords=Integer.parseInt(b);
			assert(maxRecords>=1) : "Max records must be at least 1.";
		}else if(a.equals("recordsperlevel")){
			recordsPerLevel=Integer.parseInt(b);
		}else if(a.equals("format")){
			assert(b!=null) : "Invalid format: "+arg;
			if(b.equalsIgnoreCase("json")){
				format=FORMAT_JSON;
			}else if(b.equalsIgnoreCase("jsonarray")){
				format=FORMAT_JSON;
				jsonArray=true;
			}else if(b.equalsIgnoreCase("d3")){
				format=FORMAT_JSON;
				printD3=true;
			}else if(b.equalsIgnoreCase("constellation")){
				format=FORMAT_CONSTELLATION;
			}else if(b.equalsIgnoreCase("3column") || b.equalsIgnoreCase("queryrefani")){
				format=FORMAT_QUERY_REF_ANI;
			}else if(Tools.isDigit(b.charAt(0))){
				format=Integer.parseInt(b);
			}else{
				assert(false) : "Invalid format: "+arg;
			}
		}else if(a.equalsIgnoreCase("json")){
			if(Parse.parseBoolean(b)){
				format=FORMAT_JSON;
			}else{
				if(format==FORMAT_JSON){format=default_format;}
			}
		}else if(a.equalsIgnoreCase("jsonarray") || a.equalsIgnoreCase("jsonarrays")){
			if(Parse.parseBoolean(b)){
				format=FORMAT_JSON;
				jsonArray=true;
			}else{
				jsonArray=false;
			}
		}else if(a.equalsIgnoreCase("d3") || a.equalsIgnoreCase("printd3")){
			if(Parse.parseBoolean(b)){
				format=FORMAT_JSON;
				printD3=true;
			}else{
				printD3=false;
			}
		}else if(a.equalsIgnoreCase("jsonarray") || a.equalsIgnoreCase("jsonarrays")){
			if(Parse.parseBoolean(b)){
				jsonArray=true;
			}else{
				jsonArray=false;
			}
		}else if(a.equalsIgnoreCase("d3levelnodes")){
			D3LevelNodes=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("d3hitsize")){
			if(Parse.parseBoolean(b)){D3sizeMode=D3_HIT_SIZE;}
		}else if(a.equalsIgnoreCase("d3anisize")){
			if(Parse.parseBoolean(b)){D3sizeMode=D3_ANI_SIZE;}
		}else if(a.equalsIgnoreCase("d3wkidsize")){
			if(Parse.parseBoolean(b)){D3sizeMode=D3_WKID_SIZE;}
		}else if(a.equalsIgnoreCase("d3depthsize")){
			if(Parse.parseBoolean(b)){
				D3sizeMode=D3_DEPTH_SIZE;
				printDepth=true;
			}
		}else if(a.equalsIgnoreCase("d3kidsize")){
			if(Parse.parseBoolean(b)){D3sizeMode=D3_KID_SIZE;}
		}else if(a.equalsIgnoreCase("D3sizeMode")){
			D3sizeMode=Integer.parseInt(b);
		}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl") || a.equals("minlevel")){
			taxLevel=TaxTree.parseLevel(b);//TODO: Change to extended
		}
		
		else if(a.equalsIgnoreCase("requireSSU")){
			requireSSU=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("minRefSizeEstimate") || a.equalsIgnoreCase("minRefSize")){
			minRefSizeEstimate=Long.parseLong(b);
		}else if(a.equalsIgnoreCase("minRefSizeBases")){
			minRefSizeBases=Long.parseLong(b);
		}
		
		else if(a.equalsIgnoreCase("printtax") || a.equalsIgnoreCase("printtaxa")){
			printTax=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printssu") || a.equalsIgnoreCase("print16s") || a.equalsIgnoreCase("ssu")){
			printSSU=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printSSULen") || a.equalsIgnoreCase("print16slen") || a.equalsIgnoreCase("ssulen")){
			printSSULen=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printssusequence") || a.equalsIgnoreCase("print16ssequence")){
			printSSUSequence=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printqueryfilename") || a.equalsIgnoreCase("printqfname") || a.equalsIgnoreCase("printqfile") || a.equalsIgnoreCase("qfname")){
			printQueryFileName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printreffilename") || a.equalsIgnoreCase("printrfname") || a.equalsIgnoreCase("printrfile") || a.equalsIgnoreCase("rfname")){
			printRefFileName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printfilename") || a.equalsIgnoreCase("printfname") || a.equalsIgnoreCase("printfile")){
			printQueryFileName=printRefFileName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printoriginalname") || a.equalsIgnoreCase("printseqname") || a.equalsIgnoreCase("printname0") || a.equals("pn0")){
			printOriginalName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printimg")){
			printImg=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printcompleteness") || a.equalsIgnoreCase("completeness") || a.equalsIgnoreCase("printcomplt")){
			printCompleteness=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printani") || a.equalsIgnoreCase("ani")){
			printAni=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printkid") || a.equalsIgnoreCase("kid")){
			printKID=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printwkid") || a.equalsIgnoreCase("wkid")){
			printWKID=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printscore") || a.equalsIgnoreCase("score")){
			printScore=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printevalue") || a.equalsIgnoreCase("evalue")){
			printEValue=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("trackcounts")){
			trackCounts=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth") || a.equalsIgnoreCase("depth")){
			printDepth=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printdepth2") || a.equalsIgnoreCase("depth2")){
			printDepth2=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("actualdepth") || a.equalsIgnoreCase("printactualdepth")){
			printActualDepth=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printvolume") || a.equalsIgnoreCase("volume")){
			printVolume=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printavgrefhits") || a.equalsIgnoreCase("printrefhits") || a.equalsIgnoreCase("avgrefhits") || a.equalsIgnoreCase("refhits")){
			printRefHits=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("sortByDepth")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.depthComparator;}
		}else if(a.equalsIgnoreCase("sortByDepth2")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.depth2Comparator;}
		}else if(a.equalsIgnoreCase("sortByVolume")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.volumeComparator;}
		}else if(a.equalsIgnoreCase("sortByScore")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.scoreComparator;}
		}
		else if(a.equalsIgnoreCase("sortByKID")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.KIDComparator;}
		}else if(a.equalsIgnoreCase("sortByWKID") || a.equalsIgnoreCase("sortByANI")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.WKIDComparator;}
		}else if(a.equalsIgnoreCase("sortBySSU") || a.equalsIgnoreCase("sortBy16S")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.SSUComparator;}
		}else if(a.equalsIgnoreCase("sortByHits") || a.equalsIgnoreCase("sortByMatches")){
			boolean x=Parse.parseBoolean(b);
			if(x){comparator=Comparison.HitsComparator;}
		}
		
		else if(a.equalsIgnoreCase("printUMatches") || a.equalsIgnoreCase("printUHits") || a.equalsIgnoreCase("printUnique")){
			printUnique=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches2") || a.equalsIgnoreCase("printUnique2") || a.equalsIgnoreCase("unique2")){
			printUnique2=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUMatches3") || a.equalsIgnoreCase("printUnique3") || a.equalsIgnoreCase("unique3")){
			printUnique3=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printUContam")){
			printUContam=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printNoHit")){
			printNoHit=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits") || a.equalsIgnoreCase("contam") || a.equalsIgnoreCase("printcontam")){
			printContam=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("contamhits2") || a.equalsIgnoreCase("contam2") || a.equalsIgnoreCase("printcontam2")){
			if(b==null || b.length()<1){
				printContam2=true;
			}else if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}else{
				printContam2=Parse.parseBoolean(b);
			}
		}else if(a.equalsIgnoreCase("contamLevel")){
			if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				contamLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				printContam2=true;
			}else if(TaxTree.levelMapExtendedContains(b)){
				contamLevel=TaxTree.stringToLevelExtended(b);
				printContam2=true;
			}
		}
		
		else if(a.equalsIgnoreCase("reportAniOnly") || a.equalsIgnoreCase("AniOnly")){
			reportAniOnly=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printMatches")){
			printMatches=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printLength")){
			printLength=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxID")){
			printTaxID=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSize")){
			printGSize=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("gSizeKMG")){
			gSizeKMG=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGC")){
			printGC=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGKmers")){
			printGKmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printCommonAncestor") || a.equalsIgnoreCase("printCA")){
			printCommonAncestor=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printCommonAncestorLevel") || a.equalsIgnoreCase("printCAL")){
			printCommonAncestorLevel=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printTaxName")){
			printTaxName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGSeqs")){
			printGSeqs=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printGBases")){
			printGBases=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("minEntropy") || a.equalsIgnoreCase("entropy") || a.equalsIgnoreCase("efilter")){
			minEntropy=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("minprob") || a.equalsIgnoreCase("pfilter")){
			minProb=(float)Double.parseDouble(b);
		}else if(a.equalsIgnoreCase("minQual") || a.equalsIgnoreCase("minq")){
			minQual=Byte.parseByte(b);
		}
		
		else if(a.equalsIgnoreCase("printColors") || a.equalsIgnoreCase("colors") || a.equalsIgnoreCase("color")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(b==null || b.length()<1){
				printColors=true;
			}else if(b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
				printColors=true;
			}else if(b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")){
				printColors=false;
			}else{
				printColors=true;
				if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
					colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
				}else{
					colorLevel=TaxTree.stringToLevelExtended(b);
				}
			}
			setColors=true;
//			System.err.println("Parsed "+arg); //123
		}else if(a.equalsIgnoreCase("colorLevel")){
//			System.err.println("Parsing '"+arg+"'"); //123
			if(Tools.isDigit(b.charAt(0)) || b.charAt(0)=='-'){
				colorLevel=Tools.max(0, TaxTree.levelToExtended(Integer.parseInt(b)));
			}else{
				colorLevel=TaxTree.stringToLevelExtended(b);
			}
//			System.err.println("Parsed "+arg); //123
		}
		
		else if(a.equalsIgnoreCase("printRefDivisor") || a.equalsIgnoreCase("printRDiv")){
			printRefDivisor=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQueryDivisor") || a.equalsIgnoreCase("printQDiv")){
			printQueryDivisor=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printRefSize") || a.equalsIgnoreCase("printRSize")){
			printRefSize=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printQuerySize") || a.equalsIgnoreCase("printQSize")){
			printQuerySize=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("printContamHits") || a.equalsIgnoreCase("printCHits")){
			printContamHits=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printIntersection") || a.equalsIgnoreCase("intersection") || a.equalsIgnoreCase("intersect")){
			printIntersection=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("mergePairs") || a.equalsIgnoreCase("merge")){
			mergePairs=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("printAll")){
			if(Parse.parseBoolean(b)){
				setPrintAll();
			}
		}
		
		else if(a.equals("samplerate")){
			samplerate=Float.parseFloat(b);
		}else if(a.equals("reads")){
			maxReads=Parse.parseKMG(b);
		}else if(a.equals("mode") || a.equalsIgnoreCase("single") || a.equalsIgnoreCase("singlesketch") || a.equalsIgnoreCase("onesketch")
				|| a.equalsIgnoreCase("persequence") || a.equalsIgnoreCase("sequence") || a.equalsIgnoreCase("pertaxa") 
				|| a.equalsIgnoreCase("perheader") || a.equalsIgnoreCase("perfile")){
			mode=SketchObject.parseMode(arg, a, b);
		}
		
		//For format 3
		else if(a.equalsIgnoreCase("useTaxidName") || a.equalsIgnoreCase("useTaxidAsName")){
			useTaxidName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useImgName") || a.equalsIgnoreCase("useImgAsName")){
			useImgName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useTaxName") || a.equalsIgnoreCase("useTaxAsName")){
			useTaxName=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("useFilePrefixName") || a.equalsIgnoreCase("useFilePrefixAsName")){
			useFilePrefixName=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("taxfilterincludelevel") || a.equalsIgnoreCase("includelevel") 
				|| a.equalsIgnoreCase("taxlevelwhite") || a.equalsIgnoreCase("ilevel") || a.equalsIgnoreCase("whitelevel")){
			taxLevelWhite=TaxTree.parseLevel(b);//TODO:  Change to extended
		}else if(a.equalsIgnoreCase("taxfilterinclude") || a.equalsIgnoreCase("include") || a.equalsIgnoreCase("taxfilterwhitelist")){
			taxFilterWhiteList=b;
		}else if(a.equalsIgnoreCase("taxfilterincludestring") || a.equalsIgnoreCase("includestring")
				|| a.equalsIgnoreCase("taxfilterwhitestring") || a.equalsIgnoreCase("istring")){
			taxFilterWhiteString=b;
		}else if(a.equalsIgnoreCase("banUnclassified") || a.equalsIgnoreCase("noUnclassified")){
			banUnclassified=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("banVirus") || a.equalsIgnoreCase("noVirus") || a.equalsIgnoreCase("banViruses") || a.equalsIgnoreCase("noViruses")){
			banVirus=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("taxfilterexcludelevel") || a.equalsIgnoreCase("excludelevel") 
				|| a.equalsIgnoreCase("taxlevelblack") || a.equalsIgnoreCase("elevel") || a.equalsIgnoreCase("blacklevel")){
			taxLevelBlack=TaxTree.parseLevel(b);//TODO:  Change to extended
		}else if(a.equalsIgnoreCase("taxfilterexclude") || a.equalsIgnoreCase("exclude") || a.equalsIgnoreCase("taxfilterblacklist")){
			taxFilterBlackList=b;
		}else if(a.equalsIgnoreCase("taxfilterexcludestring") || a.equalsIgnoreCase("excludestring")
				|| a.equalsIgnoreCase("taxfilterblackstring") || a.equalsIgnoreCase("estring")){
			taxFilterBlackString=b;
		}
		
		else if(a.equalsIgnoreCase("minkmercount") || a.equalsIgnoreCase("minkeycount") || a.equalsIgnoreCase("mincount") || a.equalsIgnoreCase("minKeyOccuranceCount")){
			minKeyOccuranceCount=Tools.max(1, Integer.parseInt(b));
		}
		
		//TODO:  Eventually remove support for "amino" and "k" and just support "hamino" and "hk"
		//This stands for "header amino" and "header k".
		
		//Parameters for compatibility verification
		else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("hk")){
//			System.err.println("A: k="+k+", k2="+k2+", arg="+arg);
			if(b.indexOf(',')>=0){
				String[] split=b.split(",");
				assert(split.length==2) : "\nBad argument "+arg+"\n"+b+"\n";
				int x=Integer.parseInt(split[0]);
				int y=Integer.parseInt(split[1]);
				k=Tools.max(x, y);
				k2=Tools.min(x, y);
				if(k==k2){k2=0;}
//				System.err.println("B: k="+k+", k2="+k2+", split="+Arrays.toString(split));
			}else{
				k=Integer.parseInt(b);
//				System.err.println("C: k="+k+", k2="+k2);
			}
		}else if(a.equalsIgnoreCase("hashversion") || a.equalsIgnoreCase("hv")){
			hashVersion=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("amino") || a.equalsIgnoreCase("hamino")){
			amino=Parse.parseBoolean(b);
			if(amino){translate=false;}
		}else if(a.equalsIgnoreCase("translate")){
			translate=Parse.parseBoolean(b);
			if(translate){amino=false;}
		}else if(a.equalsIgnoreCase("sixframes")){
			sixframes=Parse.parseBoolean(b);
			if(sixframes){amino=false; translate=true;}
		}
		
		else if(a.equalsIgnoreCase("requiredmeta") || a.equalsIgnoreCase("rmeta")){
			if(b==null){requiredMeta=null;}
			else{
				String[] split2=b.split(",");
				requiredMeta=new ArrayList<String>(split2.length);
				for(String mt : split2){
					assert(mt.indexOf(':')>=0) : "Metadata tags must contain ':' symbol: "+mt;
					requiredMeta.add(mt);
				}
			}
		}else if(a.equalsIgnoreCase("bannedmeta") || a.equalsIgnoreCase("bmeta")){
			if(b==null){bannedMeta=null;}
			else{
				String[] split2=b.split(",");
				bannedMeta=new ArrayList<String>(split2.length);
				for(String mt : split2){
					assert(mt.indexOf(':')>=0) : "Metadata tags must contain ':' symbol: "+mt;
					bannedMeta.add(mt);
				}
			}
		}
		
//		else if(a.equalsIgnoreCase("requiredtaxid") || a.equalsIgnoreCase("rtaxid")){
//			if(b==null){requiredTaxid=null;}
//			else{
//				String[] split2=b.split(",");
//				requiredTaxid=new IntList(split2.length);
//				for(String mt : split2){
//					requiredTaxid.add(Integer.parseInt(mt));
//				}
//				if(requiredTaxid.isEmpty()){requiredTaxid=null;}
//			}
//		}else if(a.equalsIgnoreCase("bannedtaxid") || a.equalsIgnoreCase("btaxid")){
//			if(b==null){bannedTaxid=null;}
//			else{
//				String[] split2=b.split(",");
//				bannedTaxid=new IntList(split2.length);
//				for(String mt : split2){
//					bannedTaxid.add(Integer.parseInt(mt));
//				}
//				if(bannedTaxid.isEmpty()){bannedTaxid=null;}
//			}
//		}
		
		else if(a.equalsIgnoreCase("requiredmetaand") || a.equalsIgnoreCase("rmetaand")){
			requiredMetaAnd=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("requiredmetaor") || a.equalsIgnoreCase("rmetaor")){
			requiredMetaAnd=!Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("bbversion")){
			inputVersion=b;
		}
		
		else{
			return false;
		}
		return true;
	}
	
	public void postParse(boolean requireTree, boolean makeTaxFilters){
		assert(!postParsed);
		synchronized(this){
			if(postParsed){return;}
			
			if(makeTaxFilters){
				if(taxFilterWhiteList!=null || taxFilterWhiteString!=null){
					taxFilterWhite=new TaxFilter(SketchObject.taxtree, true);
					taxFilterWhite.setLevel(taxLevelWhite, false);
					taxFilterWhite.makeSet();
					taxFilterWhite.addNamesOrNumbers(taxFilterWhiteList, false);
					taxFilterWhite.setContainsString(taxFilterWhiteString);
					if(requireTree){
						assert(SketchObject.taxtree!=null) : "No taxtree loaded.";
						taxFilterWhite.setTree(SketchObject.taxtree);
						taxFilterWhite.promote();
					}
				}
				
				if(taxFilterBlackList!=null || taxFilterBlackString!=null){
					taxFilterBlack=new TaxFilter(SketchObject.taxtree, false);
					taxFilterBlack.setLevel(taxLevelBlack, false);
					taxFilterBlack.makeSet();
					taxFilterBlack.addNamesOrNumbers(taxFilterBlackList, false);
					taxFilterBlack.setContainsString(taxFilterBlackString);
					if(requireTree){
						assert(SketchObject.taxtree!=null) : "No taxtree loaded.";
						taxFilterBlack.setTree(SketchObject.taxtree);
						taxFilterBlack.promote();
					}
				}
			}
			
			noFilters=(!hasMetaFilters() && !hasTaxFilters() && !requireSSU && minRefSizeEstimate<1 && minRefSizeBases<1);
			postParsed=true;
		}
	}
	
	public boolean postParsed(){return postParsed;}
	
	@Override
	public String toString(){
		return toString(-1);
	}
	
	public String toString(int chunkNum){
		StringBuilder sb=new StringBuilder();
		sb.append("##");
		sb.append("hits=").append(minHits);
		if(chunkNum>=0){sb.append(" chunk=").append(chunkNum);}
		sb.append(" wkid=").append(Tools.format("%.5f",minWKID));
		if(minANI>0){sb.append(" id=").append(Tools.format("%.5f",minANI));}
		if(minBases>0){sb.append(" minbases=").append(minBases);}
		if(minSizeRatio>0){sb.append(" minsizeratio=").append(Tools.format("%.5f",minSizeRatio));}
		sb.append(" records=").append(maxRecords);
		if(recordsPerLevel>0){sb.append(" recordsperlevel=").append(recordsPerLevel);}
		sb.append(" format=").append(format);
		sb.append(" level=").append(taxLevel);
		if(inputVersion!=null){sb.append(" bbversion=").append(inputVersion);}
		
		if(k!=SketchObject.defaultK || k2!=0 || k!=SketchObject.k || k2!=SketchObject.k2){
			assert(k>0 && k2>=0 && k2<k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			assert(SketchObject.k>0 && SketchObject.k2>=0 && SketchObject.k2<SketchObject.k) : "Bad values for k: "+k+", "+k2+", "+SketchObject.k+", "+SketchObject.k2;
			sb.append(" hk=").append(SketchObject.k).append(',').append(SketchObject.k2);
		}
		if(SketchObject.amino){sb.append(" hamino=").append(SketchObject.amino);} //TODO: This conflicts with Parser flag
		if(SketchObject.translate){sb.append(" translate=").append(SketchObject.translate);}
		if(SketchObject.sixframes){sb.append(" sixframes=").append(SketchObject.sixframes);}
		if(SketchObject.HASH_VERSION>1){sb.append(" hashversion=").append(SketchObject.HASH_VERSION);}

		if(true){sb.append(" printSSU=").append(printSSU());}
		if(requireSSU){sb.append(" requireSSU=").append(requireSSU);}
		if(minRefSizeEstimate>0){sb.append(" minRefSizeEstimate=").append(minRefSizeEstimate);}
		if(minRefSizeBases>0){sb.append(" minRefSizeBases=").append(minRefSizeBases);}
		
		if(json()){sb.append(" printSSUSequence=").append(printSSUSequence);}
		if(printSSULen){sb.append(" printSSULen=").append(printSSULen);}
		if(true || printTax!=default_printTax){sb.append(" printTax=").append(printTax);}
//		if(true || printFileName!=default_printFileName){sb.append(" printfname=").append(printFileName);}
		if(true || printQueryFileName!=default_printQueryFileName){sb.append(" printqfname=").append(printQueryFileName);}
		if(true || printRefFileName!=default_printRefFileName){sb.append(" printrfname=").append(printRefFileName);}
		if(true || printOriginalName!=default_printOriginalName){sb.append(" pn0=").append(printOriginalName);}
		if(true || printImg!=default_printImg){sb.append(" printImg=").append(printImg);}
		if(true || printAni!=default_printAni){sb.append(" printAni=").append(printAni);}
		if(!printKID){sb.append(" printKID=").append(printKID);}
		if(!printWKID){sb.append(" printWKID=").append(printWKID);}
		if(true || printCompleteness!=default_printCompleteness){sb.append(" printCompleteness=").append(printCompleteness);}

		if(true || printUnique!=default_printUnique){sb.append(" printUMatches=").append(printUnique);}
		if(true || printUnique2!=default_printUnique2){sb.append(" printUnique2=").append(printUnique2);}
		if(true || printUnique3!=default_printUnique3){sb.append(" printUnique3=").append(printUnique3);}
		if(true || printUContam!=default_printUContam){sb.append(" printUContam=").append(printUContam);}
		if(true || printNoHit!=default_printNoHit){sb.append(" printNoHit=").append(printNoHit);}
		if(true || printContam!=default_printContam){sb.append(" contam=").append(printContam);}
		if(true){sb.append(" contam2=").append(printContam2 ? TaxTree.extendedToLevel(contamLevel)+"" : "f");}

		if(true || printScore!=default_printScore){sb.append(" printScore=").append(printScore);}
		if(true || printEValue!=default_printEValue){sb.append(" printEValue=").append(printEValue);}
		
		if(true || printDepth!=default_printDepth){sb.append(" printDepth=").append(printDepth);}
		if(true || printDepth2!=default_printDepth2){sb.append(" printDepth2=").append(printDepth2);}
		if(true || printActualDepth!=default_printActualDepth){sb.append(" printActualDepth=").append(printActualDepth);}
		if(true || printVolume!=default_printVolume){sb.append(" printVolume=").append(printVolume);}
		if(true || printRefHits!=default_printRefHits){sb.append(" printRefHits=").append(printRefHits);}
		
		if(true || printMatches!=default_printMatches){sb.append(" printMatches=").append(printMatches);}
		if(true || printLength!=default_printLength){sb.append(" printLength=").append(printLength);}
		if(true || printTaxID!=default_printTaxID){sb.append(" printTaxID=").append(printTaxID);}
		if(true || printGSize!=default_printGSize){sb.append(" printGSize=").append(printGSize);}
		if(true || gSizeKMG!=default_gSizeKMG){sb.append(" gSizeKMG=").append(gSizeKMG);}
		if(true || printGC!=default_printGC){sb.append(" printGC=").append(printGC);}
		if(true || printGKmers!=default_printGKmers){sb.append(" printGKmers=").append(printGKmers);}
		
		if(printCommonAncestor){sb.append(" printCommonAncestor=").append(printCommonAncestor);}
		if(printCommonAncestorLevel){sb.append(" printCommonAncestorLevel=").append(printCommonAncestorLevel);}
		
		if(true || printTaxName!=default_printTaxName){sb.append(" printTaxName=").append(printTaxName);}
		if(true || printGSeqs!=default_printGSeqs){sb.append(" printGSeqs=").append(printGSeqs);}
		if(true || printGBases!=default_printGBases){sb.append(" printGBases=").append(printGBases);}
		if(true || minEntropy!=default_minEntropy){sb.append(" minEntropy=").append(Tools.format("%.4f", minEntropy));}
		if(true || minProb!=default_minProb){sb.append(" minProb=").append(Tools.format("%.4f", minProb));}
		if(true || minQual!=default_minQual){sb.append(" minQual=").append((int)minQual);}
		if(jsonArray!=default_jsonArray){sb.append(" jsonArray=").append(jsonArray);}
		if(printD3!=default_printD3){sb.append(" d3=").append(printD3);}
		if(printD3){
			sb.append(" D3sizeMode=").append(D3sizeMode);
			sb.append(" D3LevelNodes=").append(D3LevelNodes);
		}
		if(comparator!=Comparison.scoreComparator){sb.append(" ").append(comparator.toString());}
		
		if(taxFilterWhiteList!=null || taxFilterWhiteString!=null){
			if(taxFilterWhiteList!=null){sb.append(" taxfilterwhitelist=").append(taxFilterWhiteList);}
			if(taxFilterWhiteString!=null){sb.append(" taxfilterwhitestring=").append(taxFilterWhiteString);}
			sb.append(" taxlevelwhite=").append(taxLevelWhite);
		}
		if(taxFilterBlackList!=null || taxFilterBlackString!=null){
			if(taxFilterBlackList!=null){sb.append(" taxfilterblacklist=").append(taxFilterBlackList);}
			if(taxFilterBlackString!=null){sb.append(" taxfilterblackstring=").append(taxFilterBlackString);}
			sb.append(" taxlevelblack=").append(taxLevelBlack);
		}
		if(banUnclassified){sb.append(" banunclassified");}
		if(banVirus){sb.append(" banvirus");}
		
		if(useTaxidName){sb.append(" useTaxidName=").append(useTaxidName);}
		if(useImgName){sb.append(" useImgName=").append(useImgName);}
		if(useTaxName){sb.append(" useTaxName=").append(useTaxName);}
		
		if(true){sb.append(" colors=").append(printColors ? TaxTree.extendedToLevel(colorLevel)+"" : "f");}
		
		if(minKeyOccuranceCount!=default_minKeyOccuranceCount){sb.append(" minKeyOccuranceCount=").append(minKeyOccuranceCount);}
		
//		if(printColors && colorLevel!=default_colorLevel){sb.append(" colorLevel=").append(TaxTree.extendedToLevel(colorLevel));}
		

		if(printRefDivisor){sb.append(" printRefDivisor=").append(printRefDivisor);}
		if(printQueryDivisor){sb.append(" printQueryDivisor=").append(printQueryDivisor);}
		if(printRefSize){sb.append(" printRefSize=").append(printRefSize);}
		if(printQuerySize){sb.append(" printQuerySize=").append(printQuerySize);}
		if(printContamHits){sb.append(" printContamHits=").append(printContamHits);}
		if(printIntersection){sb.append(" printIntersection=").append(printIntersection);}
		if(mergePairs){sb.append(" mergePairs=").append(mergePairs);}
		
		if(maxReads>-1){sb.append(" reads=").append(maxReads);}
		if(mode!=default_mode){sb.append(" mode=").append(mode);}
		if(samplerate!=default_samplerate){sb.append(" samplerate=").append(Tools.format("%.4f",samplerate));}

		if(!requiredMetaAnd){sb.append(" requiredmetaand="+requiredMetaAnd);}
		if(requiredMeta!=null && !requiredMeta.isEmpty()){
			sb.append(" rmeta=");
			for(String s : requiredMeta){
				sb.append(s);
				sb.append(',');
			}
			sb.setLength(sb.length()-1);
		}
		if(bannedMeta!=null && !bannedMeta.isEmpty()){
			sb.append(" bmeta=");
			for(String s : bannedMeta){
				sb.append(s);
				sb.append(',');
			}
			sb.setLength(sb.length()-1);
		}
//		if(requiredTaxid!=null && !requiredTaxid.isEmpty()){
//			sb.append(" rtaxid=");
//			for(int i=0; i<requiredTaxid.size; i++){
//				sb.append(requiredTaxid.get(i));
//				sb.append(',');
//			}
//			sb.setLength(sb.length()-1);
//		}
//		if(bannedTaxid!=null && !bannedTaxid.isEmpty()){
//			sb.append(" btaxid=");
//			for(int i=0; i<bannedTaxid.size; i++){
//				sb.append(bannedTaxid.get(i));
//				sb.append(',');
//			}
//			sb.setLength(sb.length()-1);
//		}
		
		sb.append('\n');
		return sb.toString();
	}
	
	public boolean compatible(){
		return SketchObject.k==k && SketchObject.k2==k2 && SketchObject.aminoOrTranslate()==aminoOrTranslate() && hashVersion==SketchObject.HASH_VERSION;
	}
	
	public void setPrintAll(){
		printSSU=true;
		printSSULen=true;
		printSSUSequence=true;
		printTax=true;
		printQueryFileName=true;
		printRefFileName=true;
		printOriginalName=true;
		printImg=true;
		printAni=true;
		printKID=true;
		printWKID=true;
		printCompleteness=true;
		printScore=true;
		printEValue=true;
		printDepth=true;
		printDepth2=true;
		printVolume=true;
		printRefHits=true;
		
		printMatches=true;
		printLength=true;
		printTaxID=true;
		printGSize=true;
		printGC=true;
		printGKmers=true;
		printTaxName=true;
		printGSeqs=true;
		printGBases=true;
		
//		printColors=true;

		printUnique=true;
		printUnique2=true;
		printUnique3=true;
		printUContam=true;
		printNoHit=true;
		printContam=true;
		printContam2=true;
		
		printRefDivisor=true;
		printQueryDivisor=true;
		printRefSize=true;
		printQuerySize=true;
		printContamHits=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             JSON             ----------------*/
	/*--------------------------------------------------------------*/
	
	public JsonObject toJson(SketchResults sr){
		JsonObject j=toJson(sr.sketch);
		if(sr.list!=null){
			int i=0;
			for(Comparison c : sr.list){
				JsonObject jc=toJson(c);
				j.add(c.name(), jc);
				i++;
				if(i>=maxRecords){break;}
			}
		}
		
		if(jsonArray){
			toJsonArrayForm(j);
		}
		
		if(printD3){
			j.add("D3", toD3(sr));
		}
		
		return j;
	}
	
	public void toJsonArrayForm(JsonObject j0){
		if(j0.jmapSize()<1){return;}
		ArrayList<Object> list1=new ArrayList<Object>(j0.jmapSize());
		Object[] keys=null;
		for(Entry<String, JsonObject> e1 : j0.jmap.entrySet()){
			JsonObject j1=e1.getValue();
			ArrayList<Object> list2=new ArrayList<Object>(j1.omapSize());
			for(Entry<String, Object> e2 : j1.omap.entrySet()){
				Object o2=e2.getValue();
				list2.add(o2);
			}
			list1.add(list2.toArray());
			if(keys==null){
				ArrayList<Object> keyList=new ArrayList<Object>(j1.omapSize());
				for(Entry<String, Object> e2 : j1.omap.entrySet()){
					Object o2=e2.getKey();
					keyList.add(o2);
				}
				keys=keyList.toArray();
			}
		}

		JsonObject title=new JsonObject();
		for(Entry<String, Object> e : j0.omap.entrySet()){
			title.add(e.getKey(), e.getValue());
		}

		j0.clearJson();
		j0.clearOmap();

		j0.add("title", title);
		j0.add("header", keys);
		j0.add("rows", list1.toArray());
	}
	
	public JsonObject toJson(Sketch sk){
		assert(format==FORMAT_JSON);
		
		JsonObject j=new JsonObject();
		j.add("Name", sk.name());
		if(dbName!=null){j.add("DB", dbName);}
		j.add("SketchLen", sk.length());
		
		j.add("Seqs", sk.genomeSequences);
		j.add("Bases", sk.genomeSizeBases);
		j.add("gSize", sk.genomeSizeEstimate());
		if(sk.baseCounts!=null){j.addLiteral("GC", sk.gc(), 3);}
		if(sk.probCorrect<1 && sk.probCorrect>0){j.add("Quality", sk.probCorrect);}
		if(sk.keyCounts!=null){
			double d=Tools.averageDouble(sk.keyCounts);
			j.add("AvgCount", d);
			j.add("Depth", Tools.observedToActualCoverage(d));
		}
		
		if(sk.imgID>0){j.add("IMG", sk.imgID);}
		if(sk.spid>0){j.add("spid", sk.spid);}
		if(sk.taxID>0 && sk.taxID<SketchObject.minFakeID){j.add("TaxID", sk.taxID);}

		if((printRefFileName) && sk.fname()!=null){j.add("file", sk.fname());}
		if(printOriginalName && sk.name0()!=null){j.add("SeqName", sk.name0());}
		
		if(sk.meta!=null){
			for(String st : sk.meta){
				int colon=st.indexOf(':');
				j.add(st.substring(0,  colon), st.substring(colon+1));
			}
		}

		if(printSSULen){
			if(sk.r16SLen()>0){j.add("16SLen", sk.r16SLen());}
			if(sk.r18SLen()>0){j.add("18SLen", sk.r18SLen());}
		}
		if(printSSUSequence){
			if(sk.r16S()!=null){j.add("16SSequence", new String(sk.r16S()));}
			if(sk.r18S()!=null){j.add("18SSequence", new String(sk.r18S()));}
		}
		
		return j;
	}
	
	public JsonObject toJson(Comparison c){
		final int tid=c.taxID;
		
		JsonObject j=new JsonObject();
		
		//Text fields
		if(printTaxName){j.add("taxName", c.taxName()==null ? "." : c.taxName());}
		
		if(printCommonAncestor){j.add("commonAncestor", c.commonAncestor());}
		if(printCommonAncestorLevel){j.add("commonAncestorLevel", c.commonAncestorLevel());}
		
		if(printRefFileName){j.add("file", c.fname()==null ? "." : c.fname());}
		if(printOriginalName){j.add("seqName", c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				j.add("taxonomy", SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical, false));
			}else{
				j.add("taxonomy", (Object)null);
			}
		}
		
		if(printWKID){j.addLiteral("WKID", 100*c.wkid(), 4);}
		if(printKID){j.addLiteral("KID", 100*c.kid(), 4);}
//		if(printSSU() && c.ssuIdentity()>0){j.addLiteral("SSU", 100*c.ssuIdentity(), 3);} //Old
		if(printSSU() && c.ssuIdentity()>0){
			j.addLiteral(c.ssuType()==18 ? "18S" : "16S", 100*c.ssuIdentity(), 3);
		}
		
		//Primary fields
		if(printAni){j.addLiteral((aminoOrTranslate() ? "AAI" : "ANI"), 100*c.ani(), 3);}
		if(printCompleteness){j.addLiteral("Complt", 100*c.completeness(), 3);}
		if(printContam){j.addLiteral("Contam", 100*c.contamFraction(), 3);}
		if(printContam2){j.addLiteral("Contam2", 100*c.contam2Fraction(), 3);}
		if(printUContam){j.addLiteral("uContam", 100*c.uContamFraction(), 3);}
		if(printScore){j.add("Score", c.score());}
		if(printEValue){j.add("E-Val", Tools.format("%5.2e", c.eValue()));}
		
		if(printDepth){j.add("Depth", c.depth(printActualDepth));}
		if(printDepth2){j.add("Depth2", c.depth2(printActualDepth));}
		if(printVolume){j.add("Volume", c.volume()+0.001);}
		if(printRefHits){j.add("RefHits", c.avgRefHits());}
		
		if(printMatches){j.add("Matches", c.hits());}
		if(printUnique){j.add("Unique", c.uHits());}
		if(printUnique2){j.add("Unique2", c.unique2());}
		if(printUnique3){j.add("Unique3", c.unique3());}
		if(printNoHit){j.add("noHit", c.noHits());}
		if(printLength){j.add("Length", c.maxDivisor());}
		if(printTaxID){j.add("TaxID", tid>=SketchObject.minFakeID ? -1 : tid);}
		if(printImg){j.add("ImgID", c.imgID());}
		if(printGBases){j.add("gBases", c.genomeSizeBases());}
		if(printGKmers){j.add("gKmers", c.genomeSizeKmers());}
		if(printGSize){j.add("gSize", c.genomeSizeEstimate());}
		if(printGSeqs){j.add("gSeqs", c.genomeSequences());}
		if(c.hasGC()){j.addLiteral("GC", c.gc(), 3);}
		
		//Raw fields
		if(printRefDivisor){j.add("rDiv", c.refDivisor());}
		if(printQueryDivisor){j.add("qDiv", c.queryDivisor());}
		if(printRefSize){j.add("rSize", c.refSize());}
		if(printQuerySize){j.add("qSize", c.querySize());}
		if(printContamHits){j.add("cHits", c.contamHits());}
		


		if(printSSULen){
			if(c.has18S()){j.add("18SLen", c.b.r18SLen());}
			/*else*/ if(c.has16S()){j.add("16SLen", c.b.r16SLen());}
		}
		if(printSSUSequence){
			if(c.has18S()){j.add("18SSequence", new String(c.b.r18S()));}
			/*else*/ if(c.has16S()){j.add("16SSequence", new String(c.b.r16S()));}
		}
		
		if(printIntersection){
			Sketch intersection=Sketch.intersection(c.a, c.b);
			j.add("intersection", intersection.toString());
		}
		
		return j;
	}

	public boolean json(){return format==FORMAT_JSON;}
	
	/*--------------------------------------------------------------*/
	/*----------------              D3              ----------------*/
	/*--------------------------------------------------------------*/
	
	public JsonObject toD3(SketchResults sr){
		if(sr==null || sr.isEmpty()){return new JsonObject("name", "no hits");}
		JsonObject root=new JsonObject("name", "life");
		root.add("level", TaxTree.LIFE_E);
		if(sr.list!=null){
			int i=0;
			for(Comparison c : sr.list){
				ArrayList<JsonObject> tax=toD3List(c);
				addToLevel(root, tax, 0);
				i++;
				if(i>=maxRecords){break;}
			}
		}
		if(D3LevelNodes){
			root=converToD3ArrayFormat_LevelNode(root);
		}else{
			root=converToD3ArrayFormat_SingleNodeRoot(root);
		}
		return root;
	}
	
	private JsonObject converToD3ArrayFormat_SingleNodeRoot(JsonObject root){
		JsonObject children=root.removeJson("children");
		if(children==null){return root;}
		Object[] array=children.toJmapArray();
		root=(JsonObject)array[0];//Life node
		
		assert(root.getString("name").equalsIgnoreCase("Life")) : root;
		return converToD3ArrayFormat_SingleNode(root);
	}
	
	private JsonObject converToD3ArrayFormat_SingleNode(JsonObject nameNode){
		Object[] levelNodes=nameNode.toJmapArray();
		if(levelNodes==null){return nameNode;}
		nameNode.clearJson();
		
		ArrayList<JsonObject> fixed=new ArrayList<JsonObject>();
		for(Object o : levelNodes){
			JsonObject levelNode=(JsonObject)o;
			String level=levelNode.getString("name");
			JsonObject children=levelNode.removeJson("children");
			if(children!=null){
				Object[] childArray=children.toJmapArray();
				for(Object o2 : childArray){
					JsonObject child=(JsonObject)o2;//Now a name node
					String name=(String)child.removeObject("name");
					child.add("name", level+": "+name);
					converToD3ArrayFormat_SingleNode(child);
					fixed.add(child);
				}
			}
		}
		Object[] children=fixed.toArray();
		nameNode.add("children", children);
		return nameNode;
	}
	
	private JsonObject converToD3ArrayFormat_LevelNode(JsonObject levelNode){
		JsonObject children=levelNode.removeJson("children");
		if(children==null){return levelNode;}
		
		Object[] array=children.toJmapArray();
		levelNode.add("children", array);
		for(Object o : array){
			converToD3ArrayFormat_NameNode((JsonObject)o);
		}
		return levelNode;
	}
	
	private JsonObject converToD3ArrayFormat_NameNode(JsonObject nameNode){
		Object[] array=nameNode.toJmapArray();
		if(array==null){return nameNode;}
		
		nameNode.clearJson();
		nameNode.add("children", array);
		for(Object o : array){
			converToD3ArrayFormat_LevelNode((JsonObject)o);
		}
		return nameNode;
	}
	
	void addToLevel(JsonObject levelNode, ArrayList<JsonObject> list, int pos){
		JsonObject jo=list.get(pos);
		int rootLevel=levelNode.getInt("level");
		int joLevel=jo.getInt("level");
		if(rootLevel==joLevel){
			assert(levelNode.getString("name").equalsIgnoreCase(jo.getString("levelname"))) : levelNode+"\n"+jo;
			addAsChild(levelNode, list, pos);
		}else{
			assert(joLevel<rootLevel) : levelNode+"\n"+jo;
			assert(false) : levelNode+"\n"+jo;
		}
	}
		
	void addAsChild(JsonObject levelNode, ArrayList<JsonObject> list, int pos){
		JsonObject children=levelNode.getJson("children");
		if(children==null){
			children=new JsonObject();
			levelNode.add("children", children);
		}
		JsonObject jo=list.get(pos);
		String taxName=jo.getString("name");
		JsonObject nameNode=children.getJson(taxName);
		if(nameNode==null){
			nameNode=new JsonObject("name", taxName);
			children.add(taxName, nameNode);
		}
		Number size=jo.getNumber("size");
		Number oldSize=nameNode.getNumber("size");
		if(size!=null && (oldSize==null || oldSize.doubleValue()<size.doubleValue())){
			nameNode.add("size", jo.getNumber("size"));
			nameNode.add("kid", jo.getNumber("kid"));
			nameNode.add("wkid", jo.getNumber("wkid"));
			nameNode.add("ani", jo.getNumber("ani"));
			nameNode.add("hits", jo.getNumber("hits"));
			nameNode.add("depth", jo.getNumber("depth"));
		}

		if(pos<list.size()-1){//recur
			jo=list.get(pos+1);
			String levelName=jo.getString("levelname");
			int level=jo.getInt("level");
			JsonObject nextLevelNode=nameNode.getJson(levelName);
			if(nextLevelNode==null){
				nextLevelNode=new JsonObject("name", levelName);
				nextLevelNode.add("level", level);
				nameNode.add(levelName, nextLevelNode);
			}
			addAsChild(nextLevelNode, list, pos+1);
		}
	}
	
	int promote(int levelE) {
		if(levelE<0){return levelE;}
		while(!TaxTree.isSimple2(levelE) && levelE<TaxTree.LIFE){
			levelE++;
		}
		return levelE;
	}
	
	public ArrayList<JsonObject> toD3List(Comparison c){
		final ArrayList<TaxNode> nodes=toTNList(c.taxID);
		ArrayList<JsonObject> list=new ArrayList<JsonObject>(nodes.size());
		for(TaxNode tn : nodes){
			JsonObject jo=new JsonObject("name", tn.name);
			int levelE=promote(tn.levelExtended);
			jo.add("level", levelE);
			jo.add("levelname", TaxTree.levelToStringExtended(levelE));
			list.add(jo);
		}
		if(list.size()>0){
			JsonObject tail=list.get(list.size()-1);
			tail.add("size", toD3Size(c));
			tail.add("kid", c.kid());
			tail.add("wkid", c.wkid());
			tail.add("ani", c.ani());
			tail.add("hits", c.hits());
			tail.add("depth", c.depth(printActualDepth));
		}
		return list;
	}
	
	private Number toD3Size(Comparison c){
		if(D3sizeMode==D3_ANI_SIZE){
			return c.ani();
		}else if(D3sizeMode==D3_KID_SIZE){
			return c.kid();
		}else if(D3sizeMode==D3_WKID_SIZE){
			return c.wkid();
		}else if(D3sizeMode==D3_HIT_SIZE){
			return c.hits();
		}else if(D3sizeMode==D3_DEPTH_SIZE){
			return c.depth(printActualDepth);
		}
		assert(false) : "Invalid D3sizeMode "+D3sizeMode;
		return c.hits();
	}
	
	public ArrayList<TaxNode> toTNList(final int tid){
		final TaxTree tree=TaxTree.getTree();
		
		final ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		int nulls=0;
		{
			TaxNode tn=tree.getNode(tid);
			if(tn.isRanked() && !tn.cellularOrganisms()){list.add(tn);}
			while(tn.pid!=tn.id){
				tn=tree.getNode(tn.pid);
				if(tn.isRanked() && !tn.cellularOrganisms()){list.add(tn);}
			}
		}
		Collections.reverse(list);
		int prevLevelE=TaxTree.LIFE;
		for(int i=0; i<list.size(); i++){
			TaxNode tn=list.get(i);
			int levelE=promote(tn.levelExtended);
			
			if(!TaxTree.isSimple2(levelE) || (i>0 && levelE>=prevLevelE)){
				list.set(i, null);
				nulls++;
			}else{prevLevelE=levelE;}
		}
		if(nulls>0){Tools.condenseStrict(list);}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Formatting          ----------------*/
	/*--------------------------------------------------------------*/

	ByteBuilder queryHeader(Sketch sk){
		ByteBuilder bb=new ByteBuilder();
		if(format>2){return bb;}
		
		String color=toColor(sk.taxID);
		if(color!=null){bb.append(color);}
		
		bb.append("\nQuery: ").append(sk.name()==null ? "." : sk.name());
		if(dbName!=null){bb.append("\tDB: ").append(dbName);}
		bb.append("\tSketchLen: ").append(sk.length());
		bb.append("\tSeqs: ").append(sk.genomeSequences).append(' ');
		bb.append("\t"+(aminoOrTranslate() ? "SeqLen" : "Bases")+": ").append(sk.genomeSizeBases);
		bb.append("\tgSize: ").append(sk.genomeSizeEstimate());
		if(sk.baseCounts!=null){bb.append("\tGC: ").append(sk.gc(), 3);}
		if(sk.probCorrect<1 && sk.probCorrect>0){bb.append("\tQuality: ").append(sk.probCorrect, 4);}
		if(sk.keyCounts!=null){
			double d=Tools.averageDouble(sk.keyCounts);
			bb.append("\tAvgCount: ").append(d, 3);
			bb.append("\tDepth: ").append(Tools.observedToActualCoverage(d), 3);
		}

		if(sk.imgID>0){bb.append("\tIMG: ").append(sk.imgID);}
		if(sk.spid>0){bb.append("\tspid: ").append(sk.spid);}
		if(sk.taxID>0 && sk.taxID<SketchObject.minFakeID){bb.append("\tTaxID: ").append(sk.taxID);}

		if(printQueryFileName && sk.fname()!=null){bb.append("\tFile: "+sk.fname());}
		if(printOriginalName && sk.name0()!=null && !sk.name0().equals(sk.name())){bb.append("\tSeqName: "+sk.name0());}
		
		if(sk.meta!=null){
			for(String st : sk.meta){
				bb.append("\t").append(st.replaceFirst(":", ": "));
			}
		}
		
		if(color!=null){bb.append(Colors.RESET);}
		
		return bb;
	}
	
	int toColorTid(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return 0;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		return tn==null || tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E) ? 0 : tn.id;
	}
	
	String toColor(final int taxID){
		if(!printColors || SketchObject.taxtree==null || taxID<=0 || taxID>=SketchObject.minFakeID){return null;}
		TaxNode tn=SketchObject.taxtree.getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<colorLevel){
			tn=SketchObject.taxtree.getNode(tn.pid);
//			System.err.println(tn);
		}
		if(tn==null){
			return null;
		}else{
			if(tn.levelExtended>=TaxTree.LIFE_E || (tn.levelExtended>colorLevel && tn.levelExtended>TaxTree.PHYLUM_E)){return Colors.WHITE;}
			else{
//				System.err.println("*"+tn.id+", "+tn.id%Colors.colorArray.length);
				return Colors.colorArray[tn.id%Colors.colorArray.length];
			}
		}
	}
	
	String header(){
		if(format==FORMAT_JSON){return null;}
		final String ani=(aminoOrTranslate() ? "AAI" : "ANI");
		if(format==FORMAT_QUERY_REF_ANI || format==FORMAT_CONSTELLATION){
			if(reportAniOnly){return "#Query\tRef\t"+ani;}
			if(format==FORMAT_QUERY_REF_ANI){
				return "#Query\tRef\t"+ani+
				"\tQSize\tRefSize\tQBases\tRBases"+
				(printTaxID ? "\tQTaxID\tRTaxID" : "")+(printKID ? "\tKID" : "")+(printWKID ? "\tWKID" : "")+
				(printSSU() ? "\tSSU" : "")+(printCommonAncestorLevel ? "\tCALevel" : "");
			}
			if(format==FORMAT_CONSTELLATION){return "#Query\tRef\tKID\tWKID\t"+ani+"\tCmplt\tQSize\tRefSize\tQBases\tRefBases";}
		}
		return columnwiseHeader();
	}
	
	String columnwiseHeader(){
		final String ani=(aminoOrTranslate() ? "AAI" : "ANI");
		
		StringBuilder sb=new StringBuilder();
		
		//Numeric fields
		if(printKID){sb.append("WKID\t");}
		if(printWKID){sb.append("KID\t");}
		if(printAni){sb.append(ani+"\t");}
		if(printSSU()){sb.append("SSU\t");}
		if(printSSULen){sb.append("SSULen\t");}
		if(printCompleteness){sb.append("Complt\t");}
		if(printContam){sb.append("Contam\t");}
		if(printContam2){sb.append("Contam2\t");}
		if(printUContam){sb.append("uContam\t");}
		if(printScore){sb.append("Score\t");}
		if(printEValue){sb.append("E-Val\t");}
		
		if(printDepth){sb.append("Depth\t");}
		if(printDepth2){sb.append("Depth2\t");}
		if(printVolume){sb.append("Volume\t");}
		if(printRefHits){sb.append("RefHits\t");}
		if(printMatches){sb.append("Matches\t");}
		if(printUnique){sb.append("Unique\t");}
		if(printUnique2){sb.append("Unique2\t");}
		if(printUnique3){sb.append("Unique3\t");}
		if(printNoHit){sb.append("noHit\t");}
		if(printLength){sb.append("Length\t");}
		if(printTaxID){sb.append("TaxID\t");}
		if(printImg){sb.append("ImgID    \t");}
		if(printGBases){sb.append("gBases\t");}
		if(printGKmers){sb.append("gKmers\t");}
		if(printGSize){sb.append("gSize\t");}
		if(printGSeqs){sb.append("gSeqs\t");}
		if(printGC){sb.append("GC\t");}
		
		
		//Raw fields
		if(printRefDivisor){sb.append("rDiv\t");}
		if(printQueryDivisor){sb.append("qDiv\t");}
		if(printRefSize){sb.append("rSize\t");}
		if(printQuerySize){sb.append("qSize\t");}
		if(printContamHits){sb.append("cHits\t");}
		
		//Text fields
		if(printCommonAncestor){sb.append("CA\t");}
		if(printCommonAncestorLevel){sb.append("CALevel\t");}
		if(printTaxName){sb.append("taxName\t");}
		if(printRefFileName){sb.append("file\t");}
		if(printOriginalName){sb.append("seqName\t");}
		if(printTax && SketchObject.taxtree!=null){sb.append("taxonomy\t");}
		
		if(sb.length()>1){sb.setLength(sb.length()-1);}//trim trailing tab
		
		return sb.toString();
	}
	
	void formatComparisonColumnwise(Comparison c, ByteBuilder bb, int prevTid){
		final int tid=c.taxID;
		boolean reset=false;
		
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){bb.append(color);}
			if(underline!=null){bb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}
		
//		sb.append(Tools.format("%.2f%%\t%.2f%%", 100*c.idMinDivisor(), 100*c.idMaxDivisor()));
		if(printWKID){bb.append(100*c.wkid(), 2).append('%').tab();}
		if(printKID){bb.append(100*c.kid(), 2).append('%');}
		
//		if(printAni){sb.append(Tools.format("\t%.2f%%", 100*c.ani()));}
//		if(printCompleteness){sb.append(Tools.format("\t%.2f%%", 100*c.completeness()));}
//		if(printContam){sb.append(Tools.format("\t%.2f%%", 100*c.contamFraction()));}
//		if(printContam2){sb.append(Tools.format("\t%.2f%%", 100*c.contam2Fraction()));}
//		if(printUContam){sb.append(Tools.format("\t%.2f%%", 100*c.uContamFraction()));}
		
		if(printAni){bb.tab().append(100*c.ani(), 2).append('%');}
		if(printSSU()){
			float id=100*c.ssuIdentity();
			if(id>0){
				bb.tab().append(id, 2).append(c.ssuType()==16 ? '%' : '*'); //This is where 16S and 18S are differentiated
			}else{
				bb.tab().append('.');
			}
		}
		if(printSSULen){
			bb.tab().append(c.ssuLen());
		}
		if(printCompleteness){bb.tab().append(100*c.completeness(), 2).append('%');}
		if(printContam){bb.tab().append(100*c.contamFraction(), 2).append('%');}
		if(printContam2){bb.tab().append(100*c.contam2Fraction(), 2).append('%');}
		if(printUContam){bb.tab().append(100*c.uContamFraction(), 2).append('%');}
		if(printScore){bb.tab().append(c.scoreS());}
		if(printEValue){bb.tab().append(Tools.format("%5.2e", c.eValue()));}
		
		if(printDepth){bb.tab().append(c.depthS(printActualDepth));}
		if(printDepth2){bb.tab().append(c.depth2S(printActualDepth));}
		if(printVolume){bb.tab().append(c.volumeS());}
		if(printRefHits){bb.tab().append(c.avgRefHitsS());}
		
		if(printMatches){bb.tab().append(c.hits());}
		if(printUnique){bb.tab().append(c.uHits());}
		if(printUnique2){bb.tab().append(c.unique2());}
		if(printUnique3){bb.tab().append(c.unique3());}
		if(printNoHit){bb.tab().append(c.noHits());}
		if(printLength){bb.tab().append( c.maxDivisor());}
		if(printTaxID){bb.tab().append(tid>=SketchObject.minFakeID ? -1 : tid);}
		if(printImg){bb.tab().append(c.imgID());}
		if(printGBases){appendKMG(c.genomeSizeBases(), bb);}
		if(printGKmers){appendKMG(c.genomeSizeKmers(), bb);}
		if(printGSize){appendKMG(c.genomeSizeEstimate(), bb);}
		if(printGSeqs){appendKMG(c.genomeSequences(), bb);}
		if(printGC){bb.tab().append(c.gc(),3);}
		
		//Raw fields
		if(printRefDivisor){bb.tab().append(c.refDivisor());}
		if(printQueryDivisor){bb.tab().append(c.queryDivisor());}
		if(printRefSize){bb.tab().append(c.refSize());}
		if(printQuerySize){bb.tab().append(c.querySize());}
		if(printContamHits){bb.tab().append(c.contamHits());}
		
		//Text fields
		if(printCommonAncestor){bb.tab().append(c.commonAncestor());}
		if(printCommonAncestorLevel){bb.tab().append(c.commonAncestorLevel());}
		if(printTaxName){bb.tab().append(c.taxName()==null ? "." : c.taxName());}
		if(printRefFileName){bb.tab().append(c.fname()==null ? "." : c.fname());}
		if(printOriginalName){bb.tab().append(c.name0()==null ? "." : c.name0());}
		if(printTax && SketchObject.taxtree!=null){
			bb.tab();
			TaxNode tn=null;
			if(tid>0 && tid<SketchObject.minFakeID){
				tn=SketchObject.taxtree.getNode(tid);
			}

			if(tn!=null){
				bb.append(SketchObject.taxtree.toSemicolon(tn, SketchObject.skipNonCanonical, false));
			}else{
				bb.append('.');
			}
		}
		if(printTaxName && !printOriginalName && !printRefFileName && c.taxName()==null && c.name0()!=null){bb.tab().append(c.name0());} //Extra column
		
		if(reset){bb.append(Colors.RESET);}
		
		bb.append('\n');
		
		if(printIntersection){
			Sketch intersection=Sketch.intersection(c.a, c.b);
			bb.append(intersection.toString());
			bb.append('\n');
		}
		
	}
	
	void appendKMG(long value, ByteBuilder bb){
		if(gSizeKMG){
			bb.tab().append(toKMG(value));
		}else{
			bb.tab().append(value);
		}
	}
	
	String toKMG(long value){
		if(value<10000000L){return Long.toString(value);}
		value+=5;
		if(value<1000000000L){return value/1000L+"K";}
		if(value<1000000000000L){return value/1000000L+"M";}
		if(value<1000000000000000L){return value/1000000000L+"G";}
		return value/1000000000000L+"T";
	}
	
	void formatComparison3Column(Comparison c, ByteBuilder sb, int prevTid){
		Sketch query=c.a;
		final long sea=Tools.max(1, c.a.genomeSizeEstimate());
		final long seb=Tools.max(1, c.b.genomeSizeEstimate());
		final long ba=Tools.max(1, c.a.genomeSizeBases);
		final long bb=Tools.max(1, c.b.genomeSizeBases);
		final String qName=format==FORMAT_CONSTELLATION ? (useFilePrefixName ? query.filePrefix() : ""+query.sketchID) : useTaxidName ? ""+query.taxID :
			useImgName ?  ""+query.imgID : useTaxName ? query.taxName() : query.name();
		final String rName=format==FORMAT_CONSTELLATION ? (useFilePrefixName ? c.b.filePrefix() : ""+c.b.sketchID) : useTaxidName ? ""+c.taxID() :
			useImgName ?  ""+c.imgID() : useTaxName ? c.taxName() : c.name();
		final int tid=c.taxID;
		boolean reset=false;
		
		sb.append(qName).append('\t');
		if(printColors){
			final int ctid=toColorTid(tid);
			final int prevCtid=toColorTid(prevTid);

			final int cnum=ctid%Colors.colorArray.length;
			final int prevCnum=prevCtid%Colors.colorArray.length;

			String color=toColor(tid);
			String underline=(printColors && cnum==prevCnum && ctid!=prevCtid && (ctid>1 && prevCtid>1) ? Colors.UNDERLINE : null);

			if(color!=null){sb.append(color);}
			if(underline!=null){sb.append(underline);}
			reset=(color!=null || underline!=null);
			
//			System.err.println((color==null ? "" : color)+(underline==null ? "" : underline)+
//					tid+", "+prevTid+";     \t"+ctid+", "+prevCtid+";     \t"+cnum+", "+prevCnum+"; \t"+((underline!=null)+"")+Colors.RESET);
//			System.err.println(color==null ? "null" : color.substring(1));
		}

//		sb.append(rName).append(Tools.format("\t%.2f\t%.3f", 100*c.ani(), sea/(float)seb));
//		sb.append(rName).append(Tools.format("\t%.2f\t%d\t%d\t%d", 100*c.ani(), sea, seb, ba));
		
		//"#Query\tRef\tKID\tWKID\tANI\tCmplt\tQSize\tRefSize\tQBases\tRefBases";

		float kid=100*c.kid();
		float wkid=100*c.wkid();
		float ani=100*c.ani();
		float complt=100*c.completeness();
		float ssu=printSSU() ? 100*c.ssuIdentity() : 0;
		
		sb.append(rName).append('\t');
		if(reportAniOnly){
			sb.append(ani, 3).append('\t');
		}else if(format==FORMAT_CONSTELLATION){
			sb.append(kid, 3).append('\t');
			sb.append(wkid, 3).append('\t');
			sb.append(ani, 3).append('\t');
			sb.append(complt, 3).append('\t');
			sb.append(sea).append('\t');
			sb.append(seb).append('\t');
//			sb.append(ba).append('\t');
//			sb.append(bb).append('\t');
		}else{
			sb.append(ani, 3).append('\t');
			sb.append(sea).append('\t');
			sb.append(seb).append('\t');
			sb.append(ba).append('\t');
			sb.append(bb).append('\t');
			if(printTaxID){sb.append(c.a.taxID).append('\t');}
			if(printTaxID){sb.append(c.b.taxID).append('\t');}
			if(printKID){sb.append(kid, 3).append('\t');}
			if(printWKID){sb.append(wkid, 3).append('\t');}
			if(printSSU()){
				if(ssu>0){
					sb.append(ssu, 3).append('\t');
				}else{
					sb.append('.').append('\t');
				}
			}
			if(printCommonAncestorLevel){sb.append(c.commonAncestorLevel()).append('\t');}
		}
		sb.setLength(sb.length()-1);
		if(reset){sb.append(Colors.RESET);}
		
		sb.append('\n');
		
//		System.err.println(sb);
	}
	
	void formatComparison(Comparison c, ByteBuilder sb, int prevTaxID){
		if(format==FORMAT_MULTICOLUMN){
			formatComparisonColumnwise(c, sb, prevTaxID);
			return;
		}else if(format==FORMAT_QUERY_REF_ANI || format==FORMAT_CONSTELLATION){
			formatComparison3Column(c, sb, prevTaxID);
			return;
		}
		String complt=(printCompleteness ? Tools.format("\tcomplt %.2f%%%%", 100*c.completeness()) : "");
		String contam=(printContam ? Tools.format("\tcontam %.2f%%%%", 100*c.contamFraction()) : "");
//		String score=(printScore ? Tools.format("\tscore %.2f", c.score2()) : "");
		String score=(printScore ? "\tscore "+c.scoreS() : "");
		String depth=(printDepth ? "\tdepth "+c.depthS(printActualDepth) : "");
		String depth2=(printDepth2 ? "\tdepth2 "+c.depth2S(printActualDepth) : "");
		String volume=(printVolume ? "\tvolume "+c.volumeS() : "");
		String ccs=complt+contam+score;
		
		if(format==FORMAT_OLD){
			sb.append(Tools.format("WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d",
					100*c.wkid(), 100*c.kid(), c.hits(), c.minDivisor())+"\ttaxID "+c.taxID()+
					(printImg ? "\timgID "+c.imgID() : "")+"\tgKmers "+c.genomeSizeKmers()+"\t"+
					(c.taxName()==null ? "." : c.taxName())+
					((printOriginalName || (c.taxName()==null && c.name0()!=null)) ? "\t"+(c.name0()==null ? "." : c.name0()) : "")+"\n");
			if(printTax && SketchObject.taxtree!=null){
				if(c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
					TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
					if(tn!=null){
						PrintTaxonomy.printTaxonomy(tn, sb, SketchObject.taxtree, TaxTree.DOMAIN, SketchObject.skipNonCanonical);
					}
				}
				sb.append('\n');
			}
		}else{
			ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
			if(SketchObject.taxtree!=null && c.taxID()>=0 && c.taxID()<SketchObject.minFakeID){
				TaxNode tn=SketchObject.taxtree.getNode(c.taxID());
				while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
					tnl.add(tn);
					tn=SketchObject.taxtree.getNode(tn.pid);
				}
			}
			
			sb.append(Tools.format("WKID %.2f%%\tKID %.2f%%"+ccs+"\tmatches %d\tcompared %d\t",
					100*c.wkid(), 100*c.kid(), c.hits(), c.minDivisor()));
			sb.append("\ttaxID ").append(c.taxID()).append('\t');
			if(printImg){sb.append("\timgID ").append(c.imgID()).append('\t');}
			sb.append(c.taxName()).append('\t');
			if(printRefFileName){sb.append(c.fname()).append('\t');}
			if(printOriginalName || (c.taxName()==null && c.name0()!=null && !printRefFileName)){sb.append(c.name0()).append('\t');}
			
			if(printTax){
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
			}
			sb.append('\n');
			
			tnl.clear();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Filtering          ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean passesFilter(Sketch sk){
		assert(postParsed);
		if(noFilters){return true;}
		return passesSSUFilter(sk) && passesSizeFilter(sk) && passesTaxFilter(sk) && passesMetaFilter(sk);
	}
	
	private boolean passesTaxFilter(Sketch sk){
		if(taxFilterWhite==null && taxFilterBlack==null){return true;}
		int id=sk.taxID;
		if(id>0){
			if(banUnclassified && SketchObject.taxtree.isUnclassified(id)){return false;}
			if(banVirus && SketchObject.taxtree.isVirus(id)){return false;}
		}
		String s=sk.name();
		return passesTaxFilter(taxFilterWhite, id, s) && passesTaxFilter(taxFilterBlack, id, s);
	}
	
	private boolean passesTaxFilter(TaxFilter filter, int id, String s){
		if(filter==null){return true;}
		if(id>0 && !filter.passesFilter(id)){return false;}
//		if(id>0 && !filter.passesFilterFast(id)){return false;}
		if(s!=null && !filter.passesFilterByNameOnly(s)){return false;}
		return true;
	}
	
	private boolean passesMetaFilter(Sketch sk){
		if(requiredMeta==null && bannedMeta==null){return true;}
		return sk.passesMeta(requiredMeta, bannedMeta, requiredMetaAnd);
	}
	
	private boolean passesSSUFilter(Sketch sk){
		return !requireSSU || sk.hasSSU();
	}
	
	private boolean passesSizeFilter(Sketch sk){
		if(minRefSizeEstimate>0 && sk.genomeSizeEstimate()<minRefSizeEstimate){return false;}
		return sk.genomeSizeBases>=minRefSizeBases;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//These are shared with SketchObject
	//They do not affect anything and are just for the server to validate remote settings.
	private int hashVersion=SketchObject.HASH_VERSION;
	private int k=SketchObject.k;
	private int k2=SketchObject.k2;
	boolean amino=SketchObject.amino;
	boolean translate=SketchObject.translate;
	boolean sixframes=SketchObject.sixframes;
	private boolean aminoOrTranslate(){return amino | translate;}
	
	boolean noFilters=false;
	boolean postParsed=false;
	
	boolean amino(){return amino;}
	
	//These are unique
	public int maxRecords=default_maxRecords;
	public int recordsPerLevel=0;
	public float minANI=0;
	public int minBases=0;
	public float minSizeRatio=0;
	public float minWKID=default_minWKID;
	public int format=default_format;
	
	/** For tracking unique SendSketch queries */
	public int chunkNum=-1;
	public int minHits=default_minHits;
	public int taxLevel=default_taxLevel;
	public int mode=default_mode;
	public float samplerate=default_samplerate;
	public long maxReads=default_maxReads;
	public int minKeyOccuranceCount=default_minKeyOccuranceCount;
	public String inputVersion=null;
	
	public String dbName=null;

	boolean hasMetaFilters(){return requiredMeta!=null || bannedMeta!=null/* || requiredTaxid!=null || bannedTaxid!=null*/;}
	boolean hasTaxFilters(){return taxFilterWhite!=null || taxFilterBlack!=null || banUnclassified || banVirus;}
	boolean requireSSU=false;
	long minRefSizeEstimate=-1;
	long minRefSizeBases=-1;
	
	boolean requiredMetaAnd=true;
	ArrayList<String> requiredMeta=null;
	ArrayList<String> bannedMeta=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Print Columns        ----------------*/
	/*--------------------------------------------------------------*/

	public boolean printKID=true;
	public boolean printWKID=true;
	public boolean printSSU=true;
	public boolean printSSULen=false;
	public boolean printSSU(){return SketchObject.processSSU && printSSU;}
	public boolean printSSUSequence=default_printSSUSequence;
	
	//For format 2
	public boolean printTax=default_printTax;
	public boolean printOriginalName=default_printOriginalName;
	public boolean printQueryFileName=default_printQueryFileName;
	public boolean printRefFileName=default_printRefFileName;
	public boolean printImg=default_printImg;
	public boolean printAni=default_printAni;
	public boolean printCompleteness=default_printCompleteness;
	public boolean printScore=default_printScore;
	public boolean printEValue=default_printEValue;

	private boolean trackCounts=default_trackCounts;
	public boolean printDepth=default_printDepth;
	public boolean printDepth2=default_printDepth2;
	public boolean printActualDepth=default_printActualDepth;
	public boolean printVolume=default_printVolume;
	public boolean printRefHits=default_printRefHits;
	
	public boolean printLength=default_printLength;
	public boolean printTaxID=default_printTaxID;
	public boolean printGSize=default_printGSize;
	public boolean printGC=default_printGC;
	public boolean gSizeKMG=default_gSizeKMG;
	public boolean printGKmers=default_printGKmers;
	public boolean printCommonAncestor=default_printCommonAncestor;
	public boolean printCommonAncestorLevel=default_printCommonAncestorLevel;
	public boolean printTaxName=default_printTaxName;
	public boolean printGSeqs=default_printGSeqs;
	public boolean printGBases=default_printGBases;
	
	public boolean jsonArray=default_jsonArray;
	public boolean printD3=default_printD3;
	public boolean D3LevelNodes=false;
	public int D3sizeMode=D3_HIT_SIZE;
	public static final int D3_HIT_SIZE=0, D3_ANI_SIZE=1, D3_KID_SIZE=2, D3_WKID_SIZE=3, D3_DEPTH_SIZE=4;
	
	public float minEntropy=default_minEntropy;
	
	//For k=32:
	//0.000095f is >=Q6 (75%); 0.0008 is >=Q7 (80%); 0.0039 is >=Q8 (84%).
	//0.002f is >=Q7.53 (82.3%)
	//0.0017f is >=Q7.44 (82.0%)
	//0.6f works better for Illumina reads but this is more robust for PacBio. 
	public float minProb=0.0008f;
	public byte minQual=0;

	public boolean printUnique=default_printUnique;
	public boolean printUnique2=default_printUnique2;
	public boolean printUnique3=default_printUnique3;
	public boolean printUContam=default_printUContam;
	public boolean printNoHit=default_printNoHit;

	public boolean printColors=default_printColors;
	public boolean setColors=false;
	public int colorLevel=default_colorLevel;
	
	/** TODO: Note this is conflated between printing %contam and calculating things based on contam hits. */
	public boolean printContam=default_printContam;
	public boolean printContam2=default_printContam2;
	private int contamLevel=default_contamLevel;
	
	/** Raw fields */
	public boolean printMatches=default_printMatches;
	
	public boolean printRefDivisor=false;
	public boolean printQueryDivisor=false;
	public boolean printRefSize=false;
	public boolean printQuerySize=false;
	public boolean printContamHits=false;
	
	public boolean mergePairs=false;
	public boolean printIntersection=false;
	
	//For format 3 or 5
	public boolean useTaxidName=false;
	public boolean useImgName=false;
	public boolean useTaxName=false;
	public boolean useFilePrefixName=false;
	public boolean reportAniOnly=false;

	public int taxLevelWhite=0;
	public int taxLevelBlack=0;

	public String taxFilterWhiteList=null;
	public String taxFilterBlackList=null;

	public String taxFilterWhiteString=null;
	public String taxFilterBlackString=null;
	
	public TaxFilter taxFilterWhite=null;
	public TaxFilter taxFilterBlack=null;

	public boolean banUnclassified=false;
	public boolean banVirus=false;

	/** Make sure the settings are consistent, for CompareSketch.
	 * This is not yet complete. */
	public boolean checkValid(){
		if(printUnique2 || printUnique3){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		if(printContam2){
			assert(contamLevel()>=TaxTree.SUBSPECIES_E);
			assert(needContamCounts());
			assert(SketchObject.makeIndex);
			assert(SketchObject.taxtree!=null);
		}
		return true;
	}
	
	public boolean trackCounts() {
		return trackCounts || printDepth || printDepth2 || printVolume 
				|| comparator!=Comparison.scoreComparator || printD3; //|| minKeyOccuranceCount>1;
	}
	
	public boolean needContamCounts() {
		return printContam || printContam2 || printContamHits || printUnique || printUnique2 || printUnique3 || printUContam || printNoHit; // || true
	}
	
	public boolean needIndex(){
		return printContam2 || printUnique2 || printUnique3;
	}

	public int contamLevel() {
		return needIndex() ? contamLevel : -1;
	}
	
	public int compare(Comparison a, Comparison b){
		return comparator.compare(a, b);
	}
	
	public Comparator<Comparison> comparator=Comparison.scoreComparator;
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int FORMAT_OLD=0, FORMAT_MULTICOLUMN=2, FORMAT_QUERY_REF_ANI=3, FORMAT_JSON=4, FORMAT_CONSTELLATION=5;
	public static final boolean default_printD3=false;
	public static final boolean default_jsonArray=false;
	
	public static final int default_maxRecords=20;
	public static final float default_minWKID=0.0001f;
	public static final int default_format=FORMAT_MULTICOLUMN;
	public static final boolean default_printSSUSequence=false;
	public static final boolean default_printTax=false;
	public static final boolean default_printOriginalName=false;
	public static final boolean default_printQueryFileName=true;
	public static final boolean default_printRefFileName=false;
	public static final boolean default_printImg=false;
	public static final boolean default_printAni=true;
	public static final boolean default_printCompleteness=true;
	public static final boolean default_printScore=false;
	public static final boolean default_printEValue=false;
	
	public static final boolean default_trackCounts=false;
	public static final boolean default_printDepth=false;
	public static final boolean default_printDepth2=false;
	public static final boolean default_printActualDepth=true;
	public static final boolean default_printVolume=false;
	public static final boolean default_printRefHits=false;

	public static final boolean default_printContam=true;
	public static final boolean default_printContam2=false;
	
	public static final boolean default_printMatches=true;
	public static final boolean default_printLength=false;
	public static final boolean default_printTaxID=true;
	public static final boolean default_printGSize=true;
	public static final boolean default_printGC=false;
	public static final boolean default_gSizeKMG=true;
	public static final boolean default_printGKmers=false;
	public static final boolean default_printCommonAncestor=false;
	public static final boolean default_printCommonAncestorLevel=false;
	public static final boolean default_printTaxName=true;
	public static final boolean default_printGSeqs=true;
	public static final boolean default_printGBases=false;

	public static final float default_minEntropy=0.66f;
	public static final float default_minEntropy_amino=0.70f;
	public static final float default_minProb=0.0008f;
	public static final byte default_minQual=0;

	public static final boolean default_printUnique=true;
	public static final boolean default_printUnique2=false;
	public static final boolean default_printUnique3=false;
	public static final boolean default_printUContam=false;
	public static final boolean default_printNoHit=false;

	public static final boolean default_printColors=true;
	public static final int default_colorLevel=TaxTree.FAMILY_E;

	public static final int default_taxLevel=TaxTree.SPECIES;
	public static final int default_contamLevel=TaxTree.GENUS_E;
	
	public static final int default_mode=SketchObject.ONE_SKETCH;
	
	public static final int default_minHits=3;
	public static final float default_samplerate=1;
	public static final long default_maxReads=-1;
	public static final int default_minKeyOccuranceCount=1;
	
}
