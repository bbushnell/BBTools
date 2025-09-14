package align2;

public final class PackedHeap {
	
	public PackedHeap(int maxSize){
		
		int len=maxSize+1;
		if((len&1)==1){len++;} //Array size is always even.
		
		CAPACITY=maxSize;
		array=new long[len];
//		queue=new PriorityQueue<T>(maxSize);
	}
	
	public boolean add(long t){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
//		queue.add(t);
		assert(size==0 || array[size]!=-1L);
		size++;
		array[size]=t;
		percDown(size);
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return true;
	}
	
	public long peek(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return -1L;}
//		assert(array[1]==queue.peek()) : size+", "+queue.size()+"\n"+
//			array[1]+"\n"+
//			array[2]+" , "+array[3]+"\n"+
//			array[4]+" , "+array[5]+" , "+array[6]+" , "+array[7]+"\n"+
//			queue.peek()+"\n";
		//assert(testForDuplicates());
		return array[1];
	}
	
	public long poll(){
		//assert(testForDuplicates());
//		assert(queue.size()==size);
		if(size==0){return -1L;}
		long t=array[1];
//		assert(t==queue.poll());
		array[1]=array[size];
		array[size]=-1L;
		size--;
		if(size>0){percUp(1);}
//		assert(queue.size()==size);
//		assert(queue.peek()==peek());
		//assert(testForDuplicates());
		return t;
	}
	
	private void percDown(int loc){
		//assert(testForDuplicates());
		assert(loc>0);
		if(loc==1){return;}

		int next=loc/2;
		final long a=array[loc];
		long b=array[next];
		
//		while(loc>1 && (a.site<b.site || (a.site==b.site && a.column<b.column))){
		while(loc>1 && a<b){
			array[loc]=b;
			loc=next;
			next=next/2;
			b=array[next];
		}
			
		array[loc]=a;
	}
	
	private void percUp(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		int next1=loc*2;
		int next2=next1+1;
		if(next1>size){return;}
		long a=array[loc];
		long b=array[next1];
		long c=array[next2];
		assert(a!=b);
		assert(b!=c);
		assert(b!=-1L);
		//assert(testForDuplicates());
		if(c==-1L || b<=c){
			if(a>b){
//			if((a.site>b.site || (a.site==b.site && a.column>b.column))){
				array[next1]=a;
				array[loc]=b;
				//assert(testForDuplicates());
				percUp(next1);
			}
		}else{
			if(a>c){
//			if((a.site>c.site || (a.site==c.site && a.column>c.column))){
				array[next2]=a;
				array[loc]=c;
				//assert(testForDuplicates());
				percUp(next2);
			}
		}
	}
	
	private void percUpIter(int loc){
		//assert(testForDuplicates());
		assert(loc>0 && loc<=size) : loc+", "+size;
		final long a=array[loc];
		//assert(testForDuplicates());

		int next1=loc*2;
		int next2=next1+1;
		
		while(next1<=size){
			
			long b=array[next1];
			long c=array[next2];
			assert(a!=b);
			assert(b!=c);
			assert(b!=-1L);
			
			if(c==-1L || b<=c){
//			if(c==-1L || (b.site<c.site || (b.site==c.site && b.column<c.column))){
				if(a>b){
//				if((a.site>b.site || (a.site==b.site && a.column>b.column))){
//					array[next1]=a;
					array[loc]=b;
					loc=next1;
				}else{
					break;
				}
			}else{
				if(a>c){
//				if((a.site>c.site || (a.site==c.site && a.column>c.column))){
//					array[next2]=a;
					array[loc]=c;
					loc=next2;
				}else{
					break;
				}
			}
			next1=loc*2;
			next2=next1+1;
		}
		array[loc]=a;
	}
	
	public boolean isEmpty(){
//		assert((size==0) == queue.isEmpty());
		return size==0;
	}
	
	public void clear(){
//		queue.clear();
//		for(int i=1; i<=size; i++){array[i]=-1L;}
		size=0;
	}
	
	public int size(){
		return size;
	}
	
	public static int tier(int x){
		int leading=Integer.numberOfLeadingZeros(x);
		return 31-leading;
	}
	
	public boolean testForDuplicates(){
		for(int i=0; i<array.length; i++){
			for(int j=i+1; j<array.length; j++){
				if(array[i]!=-1L && array[i]==array[j]){return false;}
			}
		}
		return true;
	}
	
	final long[] array;
	private final int CAPACITY;
	private int size=0;
	
}
