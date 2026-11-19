package test.prot;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import parse.LineParser1;
import prot.FamilyGateCodec;
import prot.FamilyAlignmentMetric;
import prot.FamilyReference;
import prot.FamilyReferenceCodec;
import prot.ProteinSequence;
import structures.ByteBuilder;

/** Five-line byte roundtrips and malformed input guards; run on Dori.
 * @author Sayu
 */
public final class FamilyReferenceCodecTest{
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No codec-test arguments");}
		for(long seed:new long[]{0,1,-1,Long.MIN_VALUE,Long.MAX_VALUE}){
			final FamilyReferenceCodec.Header h=new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,seed,3,2,1,
				FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING);
			final byte[] text=FamilyReferenceCodec.appendHeader(h,new ByteBuilder()).toBytes();
			final byte[][] header=lines(text); check(header.length==5,"Header is one five-line record");
			final FamilyReferenceCodec.Header parsed=FamilyReferenceCodec.parseHeader(header,0,parser());
			check(parsed.seed==seed && parsed.totalFamilies==3 && parsed.liveFamilies==2 && parsed.emptyFamilies==1,"Header field roundtrip");
			check(Arrays.equals(text,FamilyReferenceCodec.appendHeader(parsed,new ByteBuilder()).toBytes()),"Header byte roundtrip");
			for(int size=0;size<5;size++){
				final byte[][] shortHeader=Arrays.copyOf(header,size);
				reject(()->FamilyReferenceCodec.parseHeader(shortHeader,0,parser()));
			}
			for(int row=1;row<5;row++){
				final byte[][] badHeader=header.clone(); badHeader[row]=ascii("wrong");
				reject(()->FamilyReferenceCodec.parseHeader(badHeader,0,parser()));
			}
			final byte[][] legacy=header.clone(); legacy[0]=replace(header[0],1,"family_reference_v3");
			reject(()->FamilyReferenceCodec.parseHeader(legacy,0,parser()));
		}
		final char[] ddl=new char[1024]; for(int i=0;i<ddl.length;i++){ddl[i]=(char)(i*63);}
		ddl[1023]=65535;
		final float[] spectrum=new float[400]; Arrays.fill(spectrum,1.0f/400); spectrum[0]=-0.0f; spectrum[1]=2.0f/400;
		final byte[] raw=ascii("aBZJXACDE");
		final FamilyReference family=new FamilyReference("family:one",7,12,1200,0,new ProteinSequence("family:one",raw),ddl,spectrum,
			0.2f,30,-0.3f,1.05f,Float.intBitsToFloat(FamilyGateCodec.DISABLED_BITS),0.14400001f,raw);
		final byte[] wire=FamilyReferenceCodec.appendFiveLines(family,new ByteBuilder()).toBytes();
		final byte[][] record=lines(wire); check(record.length==5,"Exactly five physical lines");
		final FamilyReference copy=FamilyReferenceCodec.parseFiveLines(record,0,parser());
		check(Arrays.equals(wire,FamilyReferenceCodec.appendFiveLines(copy,new ByteBuilder()).toBytes()),"Lossless byte roundtrip including ambiguous consensus");
		check(Arrays.equals(copy.ddl,ddl) && Arrays.equals(copy.spectrum,spectrum),"Raw array parity");
		check(Float.floatToRawIntBits(copy.spectrum[0])==0x80000000,"Signed-zero float bits retained");
		check(copy.validWords==0 && copy.spectrumNorm>0,"K5 words do not count dimers");
		check(Float.floatToRawIntBits(copy.cosineAccept)==Float.floatToRawIntBits(1.05f),"Finite cosine above one survives without a second margin");
		check(FamilyGateCodec.disabled(copy.identityUpper),"Identity keeps canonical disabled sentinel");
		final FamilyReference above100=new FamilyReference("upper110",0,1,raw.length,0,new ProteinSequence("upper110",raw),ddl,spectrum,
			0.2f,100,-0.3f,1.05f,FamilyGateCodec.upper(100,100),0.2f,raw);
		check(FamilyReferenceCodec.parseFiveLines(lines(FamilyReferenceCodec.appendFiveLines(above100,new ByteBuilder()).toBytes()),0,parser()).identityUpper>100,
			"A margin above100 remains finite and unreachable, without truncating the formula");
		final byte[][] belowFloor=record.clone();
		belowFloor[1]=replace(replace(record[1],1,new ByteBuilder().appendFloatA48(10f).toString()),4,new ByteBuilder().appendFloatA48(20f).toString());
		reject(()->FamilyReferenceCodec.parseFiveLines(belowFloor,0,parser()));
		final FamilyReference empty=FamilyReference.empty("empty",8);
		testScalarBounds(family);
		testCompactSpectrum(family);
		testDdlWidths(family);
		check(empty.cosineLower==0 && empty.identityLower==0 && empty.blosumLower==0 &&
			empty.cosineAccept==0 && empty.identityUpper==0 && empty.blosumFinal==0,"Empty DTO fields are finite inert zeros");
		reject(()->new FamilyReference("bad_empty",8,0,0,0,null,null,null,Float.NaN,0,0,0,0,0));
		final byte[][] emptyLines=lines(FamilyReferenceCodec.appendFiveLines(empty,new ByteBuilder()).toBytes());
		check(emptyLines.length==5 && new String(emptyLines[4],StandardCharsets.US_ASCII).equals("NA"),"Five-line empty marker");
		check(FamilyReferenceCodec.parseFiveLines(emptyLines,0,parser()).ddl==null,"Empty payloads null");
		final float[] zeros=new float[400];
		check(new FamilyReference("live",9,1,100,0,new ProteinSequence("live","XXXX"),new char[1024],zeros,0,0,-1,0.75f,FamilyGateCodec.DISABLED_UPPER,0).spectrumNorm==0,"Live zero spectrum is not empty");
		final byte[][] offset=new byte[7][]; System.arraycopy(record,0,offset,2,5);
		check(FamilyReferenceCodec.parseFiveLines(offset,2,parser()).rank==7,"Block offset supported");
		testTiming(record,emptyLines,offset);
		reject(()->FamilyReferenceCodec.parseFiveLines(record,1,parser()));
		for(String bad:new String[]{"", "0\t0", "00", "p", "4ooooo"}){
			final byte[][] changed=record.clone(); changed[3]=replace(record[3],0,bad); reject(()->FamilyReferenceCodec.parseFiveLines(changed,0,parser()));
		}
		for(float bad:new float[]{Float.NaN,Float.POSITIVE_INFINITY,-0.1f,1.1f}){
			final byte[][] changed=record.clone(); changed[4]=replace(record[4],0,new ByteBuilder().appendFloatA48(bad).toString());
			reject(()->FamilyReferenceCodec.parseFiveLines(changed,0,parser()));
		}
		for(String bad:new String[]{"9223372036854775808","18446744073709551616","-1","x",""}){
			final byte[][] changed=record.clone(); changed[0]=replace(record[0],3,bad); reject(()->FamilyReferenceCodec.parseFiveLines(changed,0,parser()));
		}
		for(int field:new int[]{0,1,2,3,4,5}){
			final byte[][] changed=record.clone(); changed[1]=replace(record[1],field,new ByteBuilder().appendFloatA48(Float.NaN).toString());
			reject(()->FamilyReferenceCodec.parseFiveLines(changed,0,parser()));
		}
		final byte[][] wrongEmpty=emptyLines.clone(); wrongEmpty[2]=ascii("X"); reject(()->FamilyReferenceCodec.parseFiveLines(wrongEmpty,0,parser()));
		reject(()->FamilyReferenceCodec.parseHeader(ascii("#family_reference\tfamily_reference_v1\t5\t1024\t0\t0\t0\t0\tblosum_floor_additive_v2\ta48_f32"),parser()));
		reject(()->FamilyReferenceCodec.parseHeader(ascii("#family_reference\tfamily_reference_v2\t5\t1024\t0\t0\t0\t0\tblosum_floor_additive_v2\ta48_f32"),parser()));
		reject(()->new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,0,0,0,0,"blosum_floor_additive_v2",FamilyReferenceCodec.SPECTRUM_ENCODING));
		reject(()->new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,0,Integer.MAX_VALUE,Integer.MAX_VALUE,1,FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING));
		System.out.println("FAMILY_REFERENCE_CODEC_TEST_PASS checks="+checks);
	}
	/** Width is declared twice in the header and must agree with the complete payload. */
	private static void testDdlWidths(final FamilyReference original){
		final String[] schemas={FamilyReferenceCodec.SCHEMA,FamilyReferenceCodec.MONOMER_SCHEMA,
			FamilyReferenceCodec.METRIC_SCHEMA,FamilyReferenceCodec.METRIC_MONOMER_SCHEMA,
			FamilyReferenceCodec.ALIGNMENT_SCHEMA,FamilyReferenceCodec.ALIGNMENT_MONOMER_SCHEMA};
		for(int s=0;s<schemas.length;s++){
			final float[] spectrum=new float[s%2==0 ? 400 : 20];spectrum[0]=1;
			final FamilyAlignmentMetric[] metrics=s<2 ? new FamilyAlignmentMetric[]{FamilyAlignmentMetric.LEGACY} :
				new FamilyAlignmentMetric[]{
					new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_RIGHT_TRIM_LEFT_ALIGNED_QUERY,70),
					new FamilyAlignmentMetric(FamilyAlignmentMetric.Mode.PIN_LEFT_TRIM_RIGHT_ALIGNED_QUERY,70)};
			for(FamilyAlignmentMetric metric:metrics){
				for(int width:new int[]{512,1024}){
					final char[] ddl=Arrays.copyOf(original.ddl,width);ddl[width-1]=65535;
					final FamilyReference f=s<4 ? new FamilyReference(original.id,0,original.members,original.residues,
						original.validWords,original.consensus,ddl,spectrum,0.2f,30,-0.3f,1.05f,90,0.2f) :
						new FamilyReference(original.id,0,original.members,original.residues,original.validWords,
							original.consensus,ddl,spectrum,0.2f,30,-0.3f,1.05f,90,0.2f,original.consensusRaw,
							0,new double[10],FamilyReference.AlignmentCalibration.ENABLED);
					final FamilyReferenceCodec.Header h=new FamilyReferenceCodec.Header(schemas[s],5,width,0,1,1,0,
						FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,0,metric);
					final byte[] headerBytes=FamilyReferenceCodec.appendHeader(h,new ByteBuilder()).toBytes();
					final byte[][] header=lines(headerBytes);
					check(Arrays.equals(header[3],ascii("#ddl\t"+width+"\tA48")),"Header DDL label declares the actual width");
					check(FamilyReferenceCodec.parseMetricBoundHeader(header,0,parser()).width==width,"Declared reference width roundtrips");
					final byte[] bytes=FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder(),h).toBytes();
					final FamilyReference parsed=FamilyReferenceCodec.parseFiveLinesAnySpectrum(lines(bytes),0,parser(),null);
					check(Arrays.equals(bytes,FamilyReferenceCodec.appendFiveLines(parsed,new ByteBuilder(),h).toBytes()),
						"Every schema retains complete 512/1024 payload bytes, including register65535");
					final byte[][] wrong=header.clone();wrong[3]=ascii("#ddl\t"+(width==512 ? 1024 : 512)+"\tA48");
					reject(()->FamilyReferenceCodec.parseMetricBoundHeader(wrong,0,parser()));
					final FamilyReferenceCodec.Header other=new FamilyReferenceCodec.Header(schemas[s],5,width==512 ? 1024 : 512,0,1,1,0,
						FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,0,metric);
					reject(()->FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder(),other));
					if(width==512){
						reject(()->FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder()));
						reject(()->FamilyReferenceCodec.parseFiveLines(lines(bytes),0,parser()));
					}else if(s<4){
						check(Arrays.equals(bytes,FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder()).toBytes()),
							"Explicit width-aware output preserves the existing 1024 legacy writer bytes");
					}
				}
			}
		}
	}
	/** Literal binary fractions test the wire contract independently of the projection/export code. */
	private static void testCompactSpectrum(FamilyReference f){
		final float[] spectrum=new float[20];
		spectrum[0]=-0.0f;
		for(int i=1; i<19; i++){spectrum[i]=1f/32;}
		spectrum[19]=14f/32;
		final FamilyReference compact=new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,
			f.consensus,f.ddl,spectrum,f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,
			f.identityUpper,f.blosumFinal,f.consensusRaw,f.scalarMask,f.scalarBounds);
		final byte[] wire=FamilyReferenceCodec.appendFiveLines(compact,new ByteBuilder()).toBytes();
		final byte[][] record=lines(wire),legacy=lines(FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder()).toBytes());
		check(record.length==5 && parser().set(record[4]).terms()==20,"Compact record has five lines and exactly20 spectrum values");
		for(int i=0; i<4; i++){check(Arrays.equals(record[i],legacy[i]),"Compact storage preserves non-spectrum line "+i);}
		final FamilyReferenceCodec.ParseTiming timing=new FamilyReferenceCodec.ParseTiming();
		final FamilyReference parsed=FamilyReferenceCodec.parseFiveLinesAnySpectrum(record,0,parser(),timing);
		check(Arrays.equals(wire,FamilyReferenceCodec.appendFiveLines(parsed,new ByteBuilder()).toBytes()),"Compact record byte roundtrip");
		check(timing.records==1,"Compact parsing contributes one successful timing record");
		for(int i=0; i<20; i++){
			check(Float.floatToRawIntBits(parsed.spectrum[i])==Float.floatToRawIntBits(spectrum[i]),"Compact float bits retained at "+i);
		}
		reject(()->FamilyReferenceCodec.parseFiveLines(record,0,parser()));
		check(Arrays.equals(FamilyReferenceCodec.appendFiveLines(
			FamilyReferenceCodec.parseFiveLinesAnySpectrum(legacy,0,parser(),null),new ByteBuilder()).toBytes(),
			FamilyReferenceCodec.appendFiveLines(f,new ByteBuilder()).toBytes()),"Compact-aware parser preserves legacy bytes");
		final FamilyReferenceCodec.Header h=new FamilyReferenceCodec.Header(FamilyReferenceCodec.MONOMER_SCHEMA,5,1024,0,1,1,0,
			FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING);
		final byte[] headerBytes=FamilyReferenceCodec.appendHeader(h,new ByteBuilder()).toBytes();
		final byte[][] header=lines(headerBytes);
		final FamilyReferenceCodec.Header copy=FamilyReferenceCodec.parseHeader(header,0,parser());
		check(copy.spectrumWidth==20 && copy.schema.equals("family_reference_v8") && header.length==5,"Explicit compact schema and five-line header");
		check(Arrays.equals(headerBytes,FamilyReferenceCodec.appendHeader(copy,new ByteBuilder()).toBytes()),"Compact header byte roundtrip");
		final byte[][] wrongLabel=header.clone(); wrongLabel[4]=ascii("#spectrum\t400\ta48_f32");
		reject(()->FamilyReferenceCodec.parseHeader(wrongLabel,0,parser()));
		final byte[][] wrongSchema=header.clone(); wrongSchema[0]=replace(header[0],1,FamilyReferenceCodec.SCHEMA);
		reject(()->FamilyReferenceCodec.parseHeader(wrongSchema,0,parser()));
		for(int width:new int[]{0,19,21,399,401}){
			final ByteBuilder payload=new ByteBuilder();
			for(int i=0; i<width; i++){if(i>0){payload.tab();}payload.appendFloatA48(0);}
			final byte[][] malformed=record.clone(); malformed[4]=payload.toBytes();
			reject(()->FamilyReferenceCodec.parseFiveLinesAnySpectrum(malformed,0,parser(),null));
		}
		for(float bad:new float[]{-0.1f,1.1f}){
			final byte[][] malformed=record.clone(); malformed[4]=replace(record[4],0,new ByteBuilder().appendFloatA48(bad).toString());
			reject(()->FamilyReferenceCodec.parseFiveLinesAnySpectrum(malformed,0,parser(),null));
		}
		final byte[][] badSum=record.clone(); badSum[4]=replace(record[4],19,new ByteBuilder().appendFloatA48(0).toString());
		reject(()->FamilyReferenceCodec.parseFiveLinesAnySpectrum(badSum,0,parser(),null));
		final byte[][] empty=lines(FamilyReferenceCodec.appendFiveLines(FamilyReference.empty("compact_empty",0),new ByteBuilder()).toBytes());
		check(FamilyReferenceCodec.parseFiveLinesAnySpectrum(empty,0,parser(),null).spectrum==null,"Compact-aware empty records have no spectrum");
	}
	private static void testScalarBounds(FamilyReference f){
		final double[] bounds={-1.2905,5.2905,0,0,0,0,0,0,2.1094,9.8906};
		final FamilyReference gated=new FamilyReference(f.id,f.rank,30,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
			f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,f.identityUpper,f.blosumFinal,f.consensusRaw,17,bounds);
		final byte[][] wire=lines(FamilyReferenceCodec.appendFiveLines(gated,new ByteBuilder()).toBytes());
		final FamilyReference parsed=FamilyReferenceCodec.parseFiveLines(wire,0,parser());
		check(parsed.scalarMask==17 && Arrays.equals(bounds,parsed.scalarBounds),"Scalar mask and double bounds survive exact five-line roundtrip");
		for(String[] mutation:new String[][]{{"6","32"},{"7","1e309"},{"8","-2"},{"9","1"}}){
			final byte[][] bad=wire.clone(); bad[1]=replace(wire[1],Integer.parseInt(mutation[0]),mutation[1]);
			reject(()->FamilyReferenceCodec.parseFiveLines(bad,0,parser()));
		}
		final byte[][] small=wire.clone(); small[0]=replace(wire[0],3,"29");
		reject(()->FamilyReferenceCodec.parseFiveLines(small,0,parser()));
		final FamilyReferenceCodec.Header h=new FamilyReferenceCodec.Header(FamilyReferenceCodec.SCHEMA,5,1024,0,1,1,0,
			FamilyReferenceCodec.FORMULA_VERSION,FamilyReferenceCodec.SPECTRUM_ENCODING,3.8906);
		final byte[][] header=lines(FamilyReferenceCodec.appendHeader(h,new ByteBuilder()).toBytes());
		check(FamilyReferenceCodec.parseHeader(header,0,parser()).scalarZ==3.8906,"Scalar z survives header roundtrip");
		for(String value:new String[]{"-1","1e309"}){
			final byte[][] bad=header.clone(); bad[0]=replace(header[0],10,value);
			reject(()->FamilyReferenceCodec.parseHeader(bad,0,parser()));
		}
	}
	private static void testTiming(byte[][] record,byte[][] empty,byte[][] offset){
		final FamilyReferenceCodec.ParseTiming timing=new FamilyReferenceCodec.ParseTiming();
		checkTiming(timing,0,0,0,0,0,"New accumulator is zero");
		final FamilyReference ordinary=FamilyReferenceCodec.parseFiveLines(record,0,parser());
		final FamilyReference timed=FamilyReferenceCodec.parseFiveLines(record,0,parser(),timing);
		check(Arrays.equals(FamilyReferenceCodec.appendFiveLines(ordinary,new ByteBuilder()).toBytes(),
			FamilyReferenceCodec.appendFiveLines(timed,new ByteBuilder()).toBytes()),"Timing preserves every serialized field");
		check(Arrays.equals(ordinary.consensus.enc,timed.consensus.enc) && Arrays.equals(ordinary.ddl,timed.ddl) &&
			Arrays.equals(ordinary.spectrum,timed.spectrum) && Float.floatToRawIntBits(ordinary.spectrumNorm)==
			Float.floatToRawIntBits(timed.spectrumNorm),"Timing preserves decoded arrays and cached norm bits");
		final FamilyReference noTiming=FamilyReferenceCodec.parseFiveLines(record,0,parser(),null);
		check(Arrays.equals(FamilyReferenceCodec.appendFiveLines(ordinary,new ByteBuilder()).toBytes(),
			FamilyReferenceCodec.appendFiveLines(noTiming,new ByteBuilder()).toBytes()),"Explicit null timing preserves byte output");
		check(timing.records==1 && timing.metadataGatesNanos>=0 && timing.consensusDdlNanos>=0 &&
			timing.spectrumNanos>=0 && timing.dtoNormNanos>=0,"Successful live record contributes nonnegative elapsed intervals");
		final long metadata=timing.metadataGatesNanos,ddl=timing.consensusDdlNanos,
			spectrum=timing.spectrumNanos,dto=timing.dtoNormNanos;
		final FamilyReference emptyCopy=FamilyReferenceCodec.parseFiveLines(empty,0,parser(),timing);
		check(Arrays.equals(FamilyReferenceCodec.appendFiveLines(emptyCopy,new ByteBuilder()).toBytes(),
			FamilyReferenceCodec.appendFiveLines(FamilyReferenceCodec.parseFiveLines(empty,0,parser()),new ByteBuilder()).toBytes()),
			"Timed empty record preserves byte output");
		check(timing.records==2 && timing.metadataGatesNanos>=metadata && timing.dtoNormNanos>=dto &&
			timing.consensusDdlNanos==ddl && timing.spectrumNanos==spectrum,"Empty record adds no decode time");
		check(FamilyReferenceCodec.parseFiveLines(offset,2,parser(),timing).rank==ordinary.rank && timing.records==3,
			"Caller reuses accumulator across records and block offsets");
		final FamilyReferenceCodec.ParseTiming emptyOnly=new FamilyReferenceCodec.ParseTiming();
		FamilyReferenceCodec.parseFiveLines(empty,0,parser(),emptyOnly);
		check(emptyOnly.records==1 && emptyOnly.consensusDdlNanos==0 && emptyOnly.spectrumNanos==0,
			"Empty-only worker never charges decode stages");

		// Exercise failure before parsing, each decode stage, and both DTO construction paths.
		timedReject(null,0,timing); timedReject(record,1,timing);
		final byte[][] missing=record.clone(); missing[4]=null; timedReject(missing,0,timing);
		final byte[][] badMetadata=record.clone(); badMetadata[0]=replace(record[0],3,"-1"); timedReject(badMetadata,0,timing);
		final byte[][] badGate=record.clone(); badGate[1]=ascii("bad"); timedReject(badGate,0,timing);
		final byte[][] badDdl=record.clone(); badDdl[3]=replace(record[3],0,"p"); timedReject(badDdl,0,timing);
		final byte[][] badSpectrum=record.clone(); badSpectrum[4]=replace(record[4],0,"p"); timedReject(badSpectrum,0,timing);
		final byte[][] badDto=record.clone(); badDto[4]=replace(record[4],0,new ByteBuilder().appendFloatA48(Float.NaN).toString());
		timedReject(badDto,0,timing);
		final byte[][] badEmpty=empty.clone(); badEmpty[2]=ascii("X"); timedReject(badEmpty,0,timing);
		final byte[][] badEmptyDto=empty.clone(); badEmptyDto[0]=replace(empty[0],4,"1"); timedReject(badEmptyDto,0,timing);

		final FamilyReferenceCodec.ParseTiming merged=new FamilyReferenceCodec.ParseTiming();
		merged.add(timing); merged.add(emptyOnly);
		checkTiming(merged,timing.records+emptyOnly.records,timing.metadataGatesNanos+emptyOnly.metadataGatesNanos,
			timing.consensusDdlNanos+emptyOnly.consensusDdlNanos,timing.spectrumNanos+emptyOnly.spectrumNanos,
			timing.dtoNormNanos+emptyOnly.dtoNormNanos,"Merge retains exact worker sums");
		merged.add(new FamilyReferenceCodec.ParseTiming());
		check(merged.records==4 && timing.records==3 && emptyOnly.records==1,"Merge leaves source accumulators unchanged");
	}
	private static void timedReject(byte[][] record,int offset,FamilyReferenceCodec.ParseTiming timing){
		final long records=timing.records,metadata=timing.metadataGatesNanos,ddl=timing.consensusDdlNanos,
			spectrum=timing.spectrumNanos,dto=timing.dtoNormNanos;
		final IllegalArgumentException ordinary=failure(()->FamilyReferenceCodec.parseFiveLines(record,offset,parser()));
		final IllegalArgumentException timed=failure(()->FamilyReferenceCodec.parseFiveLines(record,offset,parser(),timing));
		check(ordinary.getClass()==timed.getClass() && ordinary.getMessage().equals(timed.getMessage()),
			"Timing preserves malformed-record exception and diagnostic");
		checkTiming(timing,records,metadata,ddl,spectrum,dto,"Failed record leaves accumulator unchanged");
	}
	private static IllegalArgumentException failure(Runnable run){
		try{run.run();}catch(IllegalArgumentException expected){return expected;}
		throw new AssertionError("Malformed record accepted");
	}
	private static void checkTiming(FamilyReferenceCodec.ParseTiming timing,long records,long metadata,long ddl,long spectrum,long dto,String message){
		check(timing.records==records && timing.metadataGatesNanos==metadata && timing.consensusDdlNanos==ddl &&
			timing.spectrumNanos==spectrum && timing.dtoNormNanos==dto,message);
	}
	private static byte[] replace(byte[] line,int field,String value){
		final LineParser1 lp=parser().set(line); lp.setBounds(field);
		return new ByteBuilder().append(Arrays.copyOfRange(line,0,lp.a())).append(value).append(Arrays.copyOfRange(line,lp.b(),line.length)).toBytes();
	}
	public static byte[][] lines(byte[] bytes){
		final ArrayList<byte[]> out=new ArrayList<byte[]>(); int start=0;
		for(int i=0;i<bytes.length;i++){if(bytes[i]=='\n'){out.add(Arrays.copyOfRange(bytes,start,i)); start=i+1;}}
		if(start<bytes.length){out.add(Arrays.copyOfRange(bytes,start,bytes.length));}
		return out.toArray(new byte[out.size()][]);
	}
	private static byte[] ascii(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static LineParser1 parser(){return new LineParser1('\t');}
	private static void reject(Runnable run){try{run.run();}catch(IllegalArgumentException expected){checks++;return;}throw new AssertionError("Malformed record accepted");}
	private static void check(boolean value,String message){if(!value){throw new AssertionError(message);}checks++;}
	private static int checks;
}
