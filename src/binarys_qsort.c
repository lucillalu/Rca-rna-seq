#include <stdint.h>
#include "binarys_qsort.h"
#include "read_seeding.h"

int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, int8_t k_off)
{
    int64_t l=0, r=n-1, m;
    uint32_t tmp = 0;
    range[0] = range[1] = -1;

    //printf("k_off = %d\n", k_off);

    if (k_off == 0)
    {
        while(l <= r)
        {
            m = (l + r)/2;
            if (key < v[m])
            {
                r = m - 1;
            }
            else if (key > v[m])
            {
                l = m + 1;
            }
            else
            {
                range[0] = range[1] = m;
                return 1;
            }
        }
    }
    else
    {
        while (l <= r)
        {
            m = (l+r)/2;
            tmp = v[m] >> k_off;
            if (tmp == key)
            {
                range[0] = range[1] = m;
                
                //run low bound
                int64_t sl=l, sr=m-1, sm;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[0] = sm;
                        sr = sm-1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }

                //run upper bound
                sl = m+1; sr = r;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[1] = sm;
                        sl = sm+1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }
                return 1;
            }
            else if (tmp > key) r = m - 1;
            else l = m + 1;
        }
    }

    return -1;
}


int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}

int64_t binsearch_interval_unipath(uint32_t x, uint32_t v[], uint32_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}

int compare_uniid(const void * a, const void * b)
{
    hit_seed* sm1 = (hit_seed *)a;
    hit_seed* sm2 = (hit_seed *)b;
    
    if(sm1->uid > sm2->uid)
        return 1;
    if(sm1->uid < sm2->uid)
        return -1;
    else
    {
        //according to read_pos and uni_pos_off detected if there is an inversion
        if(sm1->read_pos > sm2->read_pos)
            return 1;
        if(sm1->read_pos < sm2->read_pos)
            return -1;
        else    return 0;         
    }
}

// int compare_seedid(const void* a, const void* b)
// {
//     vertex_m* sm1 = (vertex_m *)a;
//     vertex_m* sm2 = (vertex_m *)b;

//     if(sm1->seed_id > sm2->seed_id)
//         return 1;
//     if(sm1->seed_id < sm2->seed_id)
//         return -1;
//     else
//     {
//         if(sm1->pos_n > sm2->pos_n)
//             return 1;
//         if(sm1->pos_n < sm2->pos_n)
//             return -1;
//         else    return 0;         
//     }
// }

// int compare_mem(const void *a , const void *b)
// {
//     vertex_m* vertex1 = (vertex_m *)a;
//     vertex_m* vertex2 = (vertex_m *)b;

//     if (vertex1->read_pos > vertex2->read_pos)
//         return 1;
//     else if (vertex1->read_pos < vertex2->read_pos)
//         return -1;
//     else
//         return 0;
// }

 
int compare_uniseed(const void *a , const void *b)
{
    uni_seed* seed1 = (uni_seed *)a;
    uni_seed* seed2 = (uni_seed *)b;

    if (seed1->ref_begin > seed2->ref_begin)
        return 1;
    else if (seed1->ref_begin < seed2->ref_begin)
        return -1;
    else 
    {
        if (seed1->read_begin > seed2->read_begin)
            return 1;
        else if (seed1->read_begin < seed2->read_begin)
            return -1;
        else
            return 0;
    }

}

int compare_uniseed2(const void * a, const void * b)  //can remove
{
    uni_seed* us1 = (uni_seed *)a;
    uni_seed* us2 = (uni_seed *)b;
    
    if(us1->read_begin > us2->read_begin)
        return 1;
    if(us1->read_begin < us2->read_begin)
        return -1;
    else
    {
        //according to read_pos and uni_pos_off detected if there is an inversion
        if(us1->ref_begin > us2->ref_begin)
            return 1;
        if(us1->ref_begin < us2->ref_begin)
            return -1;
        else
            return 0;         
    }
}