#include"top.hpp"
void addtest(epoint& res)
{
	epoint bucket1[65536];
	for(int j=0;j<32768;j++)
	{
#pragma HLS PIPELINE II=3
		bucket1[j]=UniPadd(bucket1[j],bucket1[j+32768],bls377_k,bls377_q,bls377_i);
	}
	res=bucket1[121];
}
void read_s(const uint256* scalars,s_stream* s_out)
{
	uint256 x;
	for(uint32 i=0;i<SCALAR_NUM;i++)
	{
#pragma HLS PIPELINE II=16
		uint256 temps=scalars[i];
		ap_uint<16> s_ary[16];
		for(int j=0;j<WINDOW_NUM;j++)
		{
			ap_uint<WINDOW_SIZE+1> temp=temps((j+1)*WINDOW_SIZE-1,j*WINDOW_SIZE);
			if(temp[WINDOW_SIZE-1]==1)
			{
				temps+=ONE<<((j+1)*WINDOW_SIZE);
			}
			s_ary[j]=temp(15,0);
		}
		for(int k=0;k<16;k++)
		{
			s_out->write(s_ary[k]);
		}
	}
}
void divide_s(s_stream* s_0,s_stream* s_1,s_stream* s_2,s_stream* s_3)
{
	for(uint32 i=0;i<SCALAR_NUM*16-1;i=i+3)
	{
#pragma HLS PIPELINE II=3
		s_1->write(s_0->read());
		s_2->write(s_0->read());
		s_3->write(s_0->read());
	}
}

void read_ep(const uint256* p_1,const uint256* p_2,const uint256* p_3,const uint256* p_4,const uint256* p_5,
		const uint256* p_6, ep_stream* p_s1, ep_stream* p_s2, ep_stream* p_s3)
{
	for(uint32 i=0;i<SCALAR_NUM*16-1;i+=3)
	{
#pragma HLS PIPELINE II=3
		epoint ptemp1;
		ptemp1.x=uint384(p_1[i])<<128+p_5[i](127,0);
		ptemp1.y=uint384(p_2[i])<<128+p_5[i](255,128);
		ptemp1.z=uint384(p_3[i])<<128+p_6[i](127,0);
		ptemp1.t=uint384(p_4[i])<<128+p_6[i](255,128);
		p_s1->write(ptemp1);
		epoint ptemp2;
		ptemp2.x=uint384(p_1[i+1])<<128+p_5[i+1](127,0);
		ptemp2.y=uint384(p_2[i+1])<<128+p_5[i+1](255,128);
		ptemp2.z=uint384(p_3[i+1])<<128+p_6[i+1](127,0);
		ptemp2.t=uint384(p_4[i+1])<<128+p_6[i+1](255,128);
		p_s2->write(ptemp2);
		epoint ptemp3;
		ptemp3.x=uint384(p_1[i+2])<<128+p_5[i+2](127,0);
		ptemp3.y=uint384(p_2[i+2])<<128+p_5[i+2](255,128);
		ptemp3.z=uint384(p_3[i+2])<<128+p_6[i+2](127,0);
		ptemp3.t=uint384(p_4[i+2])<<128+p_6[i+2](255,128);
		p_s3->write(ptemp3);
	}
}
void processor_1(ep_stream* ep_1,s_stream* s1,epoint& res1)
{
#pragma HLS ALLOCATION function instances=UniPadd limit=1
	epoint bucket1[BUCKET_NUM];
	ap_uint<2> flags1[BUCKET_NUM/2];
	ep_stream conf_p;
	s_stream conf_s;
	//init warning!!!!!!!
	for(uint32 i=0;i<S_PART;i++)
	{
#pragma HLS PIPELINE II=3
		epoint ptemp=ep_1->read();
		uint16 index=s1->read();
		if(index==0)
			continue;
		if(index[15]==1)
		{
			index=uint16(65535)-index+1;
			ptemp.y=bls377.q-ptemp.y;
		}
		uint16 indext=index-1;
		if(flags1[indext]==3)
		{
			conf_p.write(ptemp);
			conf_s.write(index);
		}
		else if(flags1[indext]==1)
		{
			flags1[indext][1]=1;
			bucket1[indext+32768]=UniPadd(ptemp,bucket1[indext+32768],bls377_k,bls377_q,bls377_i);
			flags1[indext][1]=0;
		}
		else
		{
			flags1[indext][0]=1;
			bucket1[indext]=UniPadd(ptemp,bucket1[indext],bls377_k,bls377_q,bls377_i);
			flags1[indext][0]=0;
		}
	}
	for(int ii=0;ii<10;ii++)
	{
		if(conf_s.empty())
			break;
		epoint ptemp=conf_p.read();
		uint16 index=conf_s.read();
		bucket1[index-1]=UniPadd(ptemp,bucket1[index-1],bls377_k,bls377_q,bls377_i);
	}
	for(int j=0;j<32768;j++)
	{
#pragma HLS PIPELINE II=3
		bucket1[j]=UniPadd(bucket1[j],bucket1[j+32768],bls377_k,bls377_q,bls377_i);
	}
	epoint tmp[32];
	epoint tmp1[32];
	int te=1023;
	for(int m=0;m<32;m++)
	{
		tmp1[m]=bucket1[te];
		te=te+1024;
	}
	for(uint32 k=0;k<BUCKET_NUM;k++)
	{
#pragma HLS PIPELINE II=3
		ap_uint<5> ind=k(4,0);
		ap_uint<10> high=k(15,6);
		ap_uint<16> ind1=(ind+1)*1024-1-high;
		if(k[5]==0)
		{
			tmp[ind]=UniPadd(tmp[ind],bucket1[ind1],bls377_k,bls377_q,bls377_i);
		}
		else
		{
			tmp1[ind]=UniPadd(tmp1[ind],tmp[ind],bls377_k,bls377_q,bls377_i);
		}
	}
	epoint xx;
	epoint yy;
	for(int kk=31;kk>0;kk=kk-1)
	{
#pragma HLS PIPELINE off
		xx=UniPadd(xx,tmp1[kk],bls377_k,bls377_q,bls377_i);
		yy=UniPadd(xx,yy,bls377_k,bls377_q,bls377_i);
	}
	for(int kkk=0;kkk<10;kkk++)
	{
#pragma HLS PIPELINE off
		yy=UniPadd(yy,yy,bls377_k,bls377_q,bls377_i);
	}
	for(int kkkk=0;kkkk<32;kkkk++)
	{
#pragma HLS PIPELINE off
		yy=UniPadd(yy,tmp[kkkk],bls377_k,bls377_q,bls377_i);
	}
	res1={yy.x,yy.y,yy.z,yy.t};
}
extern"C"{
void main_function(const uint256* p_1,const uint256* p_2,const uint256* p_3,const uint256* p_4,const uint256* p_5,
		const uint256* p_6,const uint256* scalars, epoint res){
#pragma HLS INTERFACE mode=m_axi bundle=gmem8 port=res
#pragma HLS INTERFACE mode=m_axi bundle=gmem7 port=scalars
#pragma HLS INTERFACE mode=m_axi bundle=gmem6 port=p_6
#pragma HLS INTERFACE mode=m_axi bundle=gmem5 port=p_5
#pragma HLS INTERFACE mode=m_axi bundle=gmem4 port=p_4
#pragma HLS INTERFACE mode=m_axi bundle=gmem3 port=p_3
#pragma HLS INTERFACE mode=m_axi bundle=gmem2 port=p_2
#pragma HLS INTERFACE mode=m_axi bundle=gmem1 port=p_1

#pragma HLS DATAFLOW
	ep_stream ep_s1,ep_s2,ep_s3;
	s_stream index_s0,index_s1,index_s2,index_s3;
	read_ep(p_1,p_2,p_3,p_4,p_5,p_6,&ep_s1,&ep_s2,&ep_s3);
	read_s(scalars,&index_s0);
	divide_s(&index_s0,&index_s1,&index_s2,&index_s3);
//	main_steps(&ep_s1,&index_s1);
//	main_steps(&ep_s2,&index_s2);
//	main_steps(&ep_s3,&index_s3);


}
}


