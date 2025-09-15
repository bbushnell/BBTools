package structures;

public abstract class AbstractBitSet {
	
	public static AbstractBitSet make(int elements, int bitsPerElement){
		assert(bitsPerElement==1 || bitsPerElement==2) : bitsPerElement;
		assert(Integer.bitCount(bitsPerElement)==1) : bitsPerElement;
		assert(Integer.bitCount(1+Integer.numberOfTrailingZeros(bitsPerElement))==1) : bitsPerElement;
		//Can also assert
		if(bitsPerElement==1){
			return new RawBitSet(elements);
		}else if(bitsPerElement==2){
			return new MultiBitSet(elements);
		}else{
			throw new RuntimeException(""+bitsPerElement);
		}
	}
	
//	public final void set(int x){increment(x);}
	public final void increment(int x){increment(x, 1);}
	public abstract void increment(int x, int incr);
	
//	public final boolean get(int x){return getCount(x)>0;}
	public abstract int getCount(int x);
	
	/** Clears the input BitSet */
	public final void add(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){add((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){add((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	/** Clears the input BitSet */
	public final void add(RawBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		RawBitSet bs2=(RawBitSet)this;
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int[] rbs2Array=bs2.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
//			if(value!=0){bs2.addToCell(i, value);}
			rbs2Array[i]|=value;
			rbsArray[i]=0;
		}
	}
	
	/** Clears the input BitSet */
	public final void add(MultiBitSet bs){
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		MultiBitSet bs2=(MultiBitSet)this;
		assert(bits()==bs.bits());
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			if(value!=0){bs2.addToCell(i, value);}
			rbsArray[i]=0;
		}
	}
	
	public final void setToMax(AbstractBitSet bs){
		if(bs.getClass()==RawBitSet.class){setToMax((RawBitSet)bs);}
		else if(bs.getClass()==MultiBitSet.class){setToMax((MultiBitSet)bs);}
		else{throw new RuntimeException("Bad class: "+bs.getClass());}
	}
	
	public void setToMax(RawBitSet bs) {
		add(bs);
	}
	
	public void setToMax(MultiBitSet bs) {
		assert(this.getClass()==bs.getClass()) : this.getClass()+", "+bs.getClass();
		assert(bits()==bs.bits());
		assert(capacity()==bs.capacity()) : capacity()+", "+bs.capacity();
		final int[] rbsArray=bs.array();
		final int rbsLength=bs.length();
		for(int i=0; i<rbsLength; i++){
			final int value=rbsArray[i];
			if(value!=0){setToMax(i, value);}
		}
	}

	public abstract void addToCell(final int cell, final int mask);
	public abstract void setToMax(final int cell, final int mask);
	
	public abstract void clear();
	public abstract void setCapacity(long capacity, int extra);
	public abstract long cardinality();
	public abstract long capacity();
	public abstract int length();
	public abstract int bits(); //per element
	
	@Override
	public final String toString(){
		
		StringBuilder sb=new StringBuilder();
		
		final long cap=capacity();
		String spacer="";
		sb.append("{");
		for(long i=0; i<cap; i++){
			int x=getCount((int)i);
			if(x>0){
				sb.append(spacer);
				sb.append("("+i+","+x+")");
				spacer=", ";
			}
		}
		sb.append("}");
		
		return sb.toString();
	}
	
//	public final RawBitSet toRaw(){
//		if(this.getClass()==RawBitSet.class){return (RawBitSet)this;}
//		final int cap=(int)capacity();
//		RawBitSet rbs=new RawBitSet(cap, 0);
//		for(int i=0; i<cap; i++){
//			if(get(i)){rbs.set(i);}
//		}
//		return rbs;
//	}
	
}
