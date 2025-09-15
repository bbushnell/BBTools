package var2;

import java.util.Arrays;

import shared.Parse;
import shared.Tools;

/**
 * Utility class for parsing VCF format data into Var objects.
 * Handles complex VCF parsing including coordinate conversion, allele normalization,
 * and extraction of statistical information from INFO fields.
 * 
 * Supports both basic and extended parsing modes:
 * - Basic: Creates Var with position and allele data only
 * - Extended: Includes full statistical data (coverage, quality, bias metrics)
 * 
 * @author Brian Bushnell
 * @author Isla Winglet
 * @date June 25, 2025
 */
public class VcfToVar {



	/**
	 * Creates variant from VCF format line with full parsing of INFO fields.
	 * Handles complex VCF parsing including allele normalization for indels,
	 * coordinate conversion, and optional parsing of coverage/extended statistics.
	 * 
	 * @param line VCF format line as byte array
	 * @param scafMap Scaffold mapping for chromosome name resolution
	 * @param parseCoverage Whether to parse coverage statistics from INFO field
	 * @param parseExtended Whether to parse extended statistical fields
	 * @return New Var object with data from VCF line
	 */
	public static Var fromVCF(byte[] line, ScafMap scafMap, boolean parseCoverage, boolean parseExtended) {
		int a=0, b=0;
		
		// Field 0: CHROM - chromosome/scaffold name
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		String scaf=new String(line, a, b-a);
		b++;
		a=b;
		
		// Field 1: POS - position (1-based in VCF)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		int pos=Parse.parseInt(line, a, b);
		b++;
		a=b;

		// Field 2: ID - variant identifier (skip)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		b++;
		a=b;
		
		// Field 3: REF - reference allele
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		int reflen=line[a]=='.' ? 0 : b-a;  // Handle missing reference
		b++;
		a=b;

		// Field 4: ALT - alternative allele
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		byte[] alt;
		if(b<=a+1){
			alt=Var.AL_MAP[line[a]];  // Single base alternative
		}else{
			alt=Arrays.copyOfRange(line, a, b);  // Multi-base alternative
		}
		b++;
		a=b;

		// Skip QUAL field, find INFO section for extended parsing
		int infoStart=b;
		
		// Convert VCF coordinates to internal format
		final int start;
		if(alt.length!=reflen && alt.length>0){  // Indel normalization
			alt=Arrays.copyOfRange(alt, 1, alt.length);  // Remove first base
			start=pos;  // Use 1-based position directly
			if(alt.length==0){alt=Var.AL_0;}  // Empty alternative for deletions
			else if(alt.length==1 && Var.AL_MAP[alt[0]]!=null){
				alt=Var.AL_MAP[alt[0]];  // Use pre-allocated single base array
			}
			if(reflen==1){reflen--;}  // Adjust reference length for insertions
		}else{
			start=pos-1;  // Convert to 0-based coordinates
		}
		
		final int readlen=(alt==null || alt.length==0 || alt[0]=='.' ? 0 : alt.length);
		final int stop=start+reflen;
		
		// Resolve scaffold name to number
		assert(scaf!=null);
		assert(scafMap!=null);
		final int scafNum=scafMap.getNumber(scaf);
		assert(scafNum>=0) : scaf+"\n"+scafMap.keySet()+"\n"+scafMap.altKeySet()+"\n";
		
		// Determine variant type
		int type=-1;
		if(parseExtended){
			type=Tools.max(type, parseVcfIntDelimited(line, "TYP=", infoStart));
			if(type<0){type=Var.typeStartStop(start, stop, alt);}
		}else{
			type=Var.typeStartStop(start, stop, alt);
		}
		
		Var v=new Var(scafNum, start, stop, alt, type);
		
		// Parse coverage statistics if requested
		if(parseCoverage){
			infoStart=Tools.indexOfDelimited(line, "R1P=", infoStart, (byte)';');
			
			// Parse strand-specific read counts: R1P=2;R1M=0;R2P=0;R2M=0;
			int r1p=Tools.max(0, parseVcfIntDelimited(line, "R1P=", infoStart));
			int r1m=Tools.max(0, parseVcfIntDelimited(line, "R1M=", infoStart));
			int r2p=Tools.max(0, parseVcfIntDelimited(line, "R2P=", infoStart));
			int r2m=Tools.max(0, parseVcfIntDelimited(line, "R2M=", infoStart));
			
			// Parse coverage: AD=2;DP=24;MCOV=0;PPC=0;
			int cov=parseVcfIntDelimited(line, "DP=", infoStart);
			assert(cov>0) : new String(line, infoStart, line.length-infoStart);
			int mcov=parseVcfIntDelimited(line, "MCOV=", infoStart);
			
			v.coverage=cov;
			v.minusCoverage=mcov;
			v.r1plus=r1p;
			v.r1minus=r1m;
			v.r2plus=r2p;
			v.r2minus=r2m;
		}
		
		// Parse extended statistical fields if requested
		if(parseExtended){
			infoStart=Tools.indexOfDelimited(line, "PPC=", infoStart, (byte)';');
			
			// Parse pairing and allele fraction: PPC=0;AF=0.0833;RAF=0.0833;LS=280;
			int pc=Tools.max(0, parseVcfIntDelimited(line, "PPC=", infoStart));
			double raf=parseVcfDoubleDelimited(line, "RAF=", infoStart);
			long ls=Tools.max(0, parseVcfLongDelimited(line, "LS=", infoStart));
	
			// Parse mapping quality: MQS=86;MQM=43;BQS=64;BQM=32;
			long mqs=Tools.max(0, parseVcfLongDelimited(line, "MQS=", infoStart));
			int mqm=Tools.max(0, parseVcfIntDelimited(line, "MQM=", infoStart));
			long bqs=Tools.max(0, parseVcfLongDelimited(line, "BQS=", infoStart));
			int bqm=Tools.max(0, parseVcfIntDelimited(line, "BQM=", infoStart));
			
			// Parse distance and identity: EDS=18;EDM=9;IDS=1984;IDM=992;
			long eds=Tools.max(0, parseVcfLongDelimited(line, "EDS=", infoStart));
			int edm=Tools.max(0, parseVcfIntDelimited(line, "EDM=", infoStart));
			long ids=Tools.max(0, parseVcfLongDelimited(line, "IDS=", infoStart));
			int idm=Tools.max(0, parseVcfIntDelimited(line, "IDM=", infoStart));
			
			// Parse annotation fields: NVC=0;FLG=0;CED=601;HMP=0;SB=0.1809;DP4=37,28,15,97
			int nvc=Tools.max(0, parseVcfIntDelimited(line, "NVC=", infoStart));
			int flg=Tools.max(0, parseVcfIntDelimited(line, "FLG=", infoStart));
			
			// Store all extended statistics in variant object
			v.properPairCount=pc;
			v.revisedAlleleFraction=raf;
			v.lengthSum=ls;
			v.mapQSum=mqs;
			v.mapQMax=mqm;
			v.baseQSum=bqs;
			v.baseQMax=bqm;
			v.endDistSum=eds;
			v.endDistMax=edm;
			v.idSum=ids;
			v.idMax=idm;
			v.nearbyVarCount=nvc;
			v.flagged=(flg>0);
		}
		
		return v;
	}
	
	/**
	 * Parses integer value from VCF INFO field with specified key.
	 * Searches for key=value pattern and extracts integer, with bounds checking.
	 * 
	 * @param line VCF line data
	 * @param query Search key (e.g., "DP=")
	 * @param start Starting position for search
	 * @return Parsed integer value, or -1 if not found
	 */
	private static int parseVcfIntDelimited(byte[] line, String query, int start){
		return (int)Tools.min(Integer.MAX_VALUE, parseVcfLongDelimited(line, query, start));
	}
	
	/**
	 * Parses long integer value from VCF INFO field with specified key.
	 * Handles negative values and validates numeric parsing.
	 * 
	 * @param line VCF line data
	 * @param query Search key (e.g., "MQS=")
	 * @param start Starting position for search
	 * @return Parsed long value, or -1 if not found
	 */
	private static long parseVcfLongDelimited(byte[] line, String query, final int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
		if(loc<0){return -1;}
		long current=0;
		long mult=1;
		if(loc>0){
			if(line[loc+query.length()]=='-'){mult=-1; loc++;}  // Handle negative values
			for(int i=loc+query.length(); i<line.length; i++){
				final byte x=line[i];
				if(Tools.isDigit(x)){
					current=current*10+(x-'0');  // Parse digit
				}else{
					assert(x==tab || x==colon) : x+", "+query+", "+new String(line, loc, i-loc+1);
					break;  // Stop at delimiter
				}
			}
		}
		return mult*current;
	}
	
	/**
	 * Parses double/float value from VCF INFO field with specified key.
	 * Extracts decimal number from key=value pattern in INFO section.
	 * 
	 * @param line VCF line data
	 * @param query Search key (e.g., "AF=")
	 * @param start Starting position for search
	 * @return Parsed double value, or -1.0 if not found
	 */
	private static double parseVcfDoubleDelimited(byte[] line, String query, int start){
		int loc=Tools.indexOfDelimited(line, query, start, (byte)';');
		if(loc<0){return -1;}
		if(loc>0){
			loc=loc+query.length();
			int loc2=loc+1;
			// Find end of numeric value
			while(loc2<line.length){
				byte b=line[loc2];
				if(b==tab || b==colon){break;}  // Stop at field separators
				loc2++;
			}
			return Parse.parseDouble(line, loc, loc2);  // Parse substring as double
		}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Tab character for field separation */
	private static final byte tab='\t';
	/** Semicolon character for INFO field separation */
	private static final byte colon=';';
}